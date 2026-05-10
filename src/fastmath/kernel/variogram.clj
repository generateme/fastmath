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
  - Unbounded models (no sill): `power`
  - Generalised families with a configurable shape order, each available as a
    two-step creator `->name` / convenience wrapper `name`: `->bessel` / `bessel`,
    `->matern` / `matern`, `->hyperspherical` / `hyperspherical`,
    `->superspherical` / `superspherical`, `->tplstable` / `tplstable`
  - RBF-based model: `rbf->variogram` wraps any suitable RBF kernel from
    `fastmath.kernel.rbf` into a semivariogram constructor

  Empirical semivariogram:

  - `empirical` — bins all observation pairs by lag distance and applies a chosen
    estimator to compute `gamma` per bin; returns a sequence of `{:h :gamma :n}` maps

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

  - `bounding-box-diagonal` — computes the diagonal of the bounding box of spatial points;
    used to derive the default cutoff in `empirical`
  - `remove-outliers-fence` — removes outliers using Tukey's fences criterion
  - `remove-outliers-mad` — removes outliers using median absolute deviation"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.distance :as dist]
            [fastmath.stats :as stats]
            [fastmath.optimization :as optim]
            [fastmath.random :as r]
            [fastmath.special :as special])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

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
      - `:beta` (double): the power exponent controlling the rate of growth; typically
        in `(0.0, 2.0)` for valid conditionally negative-definite models

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  See also [[exponential]], [[gaussian]], [[spherical]], [[expower]], [[tpower]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m/pow (m// x range) beta)
             (m/* psill)
             (m/+ nugget)))))

(defn exponential
  "Exponential semivariogram model.

  The formula is: `nugget + psill * (1 - exp(-x / range))`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the maximum contribution of the structured
        spatial variance above the nugget; must be non-negative
      - `:range` (double): the practical range parameter controlling how fast the semivariance
        rises towards the sill; must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  See also [[gaussian]], [[spherical]], [[power]], [[expower]], [[hole]]."
  [{:keys [^double nugget ^double psill ^double range]}]
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
        must be in `(0.0, 2.0]` for a valid conditionally negative-definite model

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  See also [[exponential]], [[gaussian]], [[power]], [[cauchy]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
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

  See also [[exponential]], [[expower]], [[spherical]], [[cubic]]."
  [{:keys [^double nugget ^double psill ^double range]}]
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
        must be positive for a bounded, valid model; clamped to `EPSILON` if near zero

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  See also [[cauchy]], [[gaussian]], [[expower]], [[power]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (let [b (m/copy-sign (m/max (m/abs beta) m/EPSILON) beta)
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

  See also [[linear]], [[cubic]], [[pentaspherical]], [[circular]], [[gaussian]]."
  [{:keys [^double nugget ^double psill ^double range]}]
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

  See also [[spherical]], [[pentaspherical]], [[gaussian]], [[circular]]."
  [{:keys [^double nugget ^double psill ^double range]}]
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

  See also [[spherical]], [[cubic]], [[circular]], [[linear]]."
  [{:keys [^double nugget ^double psill ^double range]}]
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

  See also [[spherical]], [[pentaspherical]], [[cubic]], [[linear]]."
  [{:keys [^double nugget ^double psill ^double range]}]
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

  See also [[spherical]], [[power]], [[tpower]], [[circular]]."
  [{:keys [^double nugget ^double psill ^double range]}]
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

  See also [[power]], [[linear]], [[spherical]], [[expower]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
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

  where `sinc(t) = sin(t) / t` (unnormalized sinc), with `sinc(0) = 1`.

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the amplitude of the oscillation around the sill;
        must be non-negative
      - `:range` (double): the scale parameter controlling the period of oscillation;
        must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  See also [[exponential]], [[cauchy]], [[bessel]]."
  [{:keys [^double nugget ^double psill ^double range]}]
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m// x range)
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

  See also [[rational]], [[gaussian]], [[exponential]], [[matern]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
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

(defn ->bessel
  "Creates a BesselJ semivariogram model constructor for a given Bessel function `order`.

  The formula is: `nugget + psill * (1 - Γ(order+1) * (2 * range / x)^order * J_order(x / range))`

  where `Γ` is the gamma function and `J_order` is the Bessel function of the first kind.

  Parameters:

  - `order` (double, optional): the order of the Bessel function of the first kind;
    must be non-negative; defaults to `0.0`

  Returns a semivariogram model constructor `(fn [params])` that accepts a map with:
  - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
    variation or measurement error; must be non-negative
  - `:psill` (double): the partial sill; the amplitude of the oscillation; must be non-negative
  - `:range` (double): the scale parameter controlling the period of oscillation; must be positive

  The constructor itself returns a function `(fn ^double [^double x])` that computes the
  semivariance for a given lag distance `x`. Returns `0.0` when `x` is zero.

  See also [[bessel]], [[->matern]], [[hole]]."
  ([] (->bessel 0.0))
  ([^double order]
   (let [g (special/gamma (m/inc order))]
     (fn [{:keys [^double nugget ^double psill ^double range]}]
       (fn ^double [^double x]
         (if (m/zero? x) 0.0
             (let [d (m// x range)]
               (->> (m/* g (m/pow (m// 2.0 d) order)
                         (special/bessel-J order d))
                    (m/- 1.0)
                    (m/* psill)
                    (m/+ nugget)))))))))

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

  See also [[->bessel]], [[matern]], [[hole]]."
  [params]
  ((->bessel (or (:order params) 1.0)) params))

(defn ->matern
  "Creates a Matérn semivariogram model constructor for a given smoothness `order`.

  The formula is: `nugget + psill * (1 - C(order) * (sqrt(2*order) * x / range)^order * K_order(sqrt(2*order) * x / range))`

  where `C(order) = 2^(1-order) / Γ(order)` is a normalisation constant ensuring the
  kernel equals `1.0` at the origin, `Γ` is the gamma function and `K_order` is the
  modified Bessel function of the second kind of the given order.

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

  See also [[matern]], [[->bessel]], [[gaussian]], [[exponential]]."
  ([] (->matern 1.0))
  ([^double order]
   (let [g (m// (m/pow 2.0 (m/- 1.0 order)) (special/gamma order))
         f (m/sqrt (m/* 2.0 order))]
     (fn [{:keys [^double nugget ^double psill ^double range]}]
       (fn ^double [^double x]
         (if (m/zero? x) 0.0
             (let [d (m/* f (m// x range))]
               (->> (m/* g (m/pow d order) (special/bessel-K order d))
                    (m/- 1.0)
                    (m/* psill)
                    (m/+ nugget)))))))))

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
        `0.5` is equivalent to the [[exponential]] model, higher values approach the [[gaussian]] model

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

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

  See also [[hyperspherical]], [[->superspherical]], [[spherical]], [[circular]]."
  ([] (->hyperspherical 2.0))
  ([^double order]
   (let [d (m/* (m/dec order) -0.5)
         denom (special/hypergeometric-2F1 0.5 d 1.5 1.0)]
     (fn [{:keys [^double nugget ^double psill ^double range]}]
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

  See also [[superspherical]], [[->hyperspherical]], [[spherical]], [[circular]], [[linear]]."
  ([] (->superspherical 1.0))
  ([^double order]
   (let [d (m/- order)
         denom (special/hypergeometric-2F1 0.5 d 1.5 1.0)]
     (fn [{:keys [^double nugget ^double psill ^double range]}]
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
    must be in `(0.0, 2.0)` for a valid conditionally negative-definite model

  The constructor itself returns a function `(fn ^double [^double x])` that computes the
  semivariance for a given lag distance `x`. Returns `0.0` when `x` is zero.

  See also [[tplstable]], [[power]], [[expower]], [[matern]]."
  ([] (->tplstable 0.5))
  ([^double H]
   (let [H2 (m/* 2.0 H)]
     (fn [{:keys [^double nugget ^double psill ^double range ^double beta]}]
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
        must be in `(0.0, 2.0)` for a valid conditionally negative-definite model
      - `:H` (double, optional): the Hurst exponent controlling long-range dependence;
        must be in `(0.0, 1.0)`; defaults to `0.5`

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  See also [[->tplstable]], [[power]], [[expower]], [[matern]]."
  [params]
  ((->tplstable (or (:H params) 0.5)) params))

(defn rbf->variogram
  "Builds a semivariogram model constructor from a Radial Basis Function (RBF) kernel.

  Converts any RBF kernel from `fastmath.kernel.rbf` into a semivariogram model by
  treating the kernel as the correlation function: semivariance equals
  `nugget + psill * (1 - kernel(x / range))`. For this transformation to produce a valid
  semivariogram, the kernel must satisfy two conditions: it must equal `1.0` at the origin
  and decay to `0.0` at infinity. The `range` parameter acts as the reciprocal of the
  RBF's `scale` shape parameter.

  The formula is: `nugget + psill * (1 - kernel(x / range))`

  Suitable kernels from `fastmath.kernel.rbf` (value `1.0` at zero, decaying to `0.0`):

  - `:gaussian` - smooth Gaussian decay
  - `:matern` - Matérn family; differentiability controlled by order
  - `:truncated-power` - compactly supported power kernel
  - `:wendland` - compactly supported, C^(2k) smooth
  - `:gneiting` - compactly supported Gneiting family
  - `:wu` - compactly supported Wu family
  - `:whittaker` - Whittaker family, normalised to `1.0` at zero
  - `:poisson` - Bessel J family; oscillates around `0.0`, use with care

  Unsuitable kernels (do not satisfy the required conditions):

  - `:thin-plate` / `:thin-plate-splines` - grows without bound
  - `:linear` - zero at origin, grows without bound
  - `:radial-powers` - zero or infinite at origin depending on exponent
  - `:generalized-multiquadratic` - not normalised to `+1.0` at zero
  - `:shifted-surface-splines` - grows with distance for positive exponent

  Parameters:

  - `kernel` (function): an RBF kernel function `(fn ^double [^double x])` returning `1.0`
    at `x = 0` and decaying to `0.0` as `x` increases; typically constructed via
    `fastmath.kernel.rbf` functions

  Returns a semivariogram model constructor `(fn [params])` that accepts a map with:
  - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
    variation or measurement error; must be non-negative
  - `:psill` (double): the partial sill; the maximum contribution of the structured
    spatial variance above the nugget; must be non-negative
  - `:range` (double): the scale parameter passed as reciprocal of the RBF `scale`;
    must be positive

  The constructor itself returns a function `(fn ^double [^double x])` that computes the
  semivariance for a given lag distance `x`. Returns `0.0` when `x` is zero.

  See also [[exponential]], [[gaussian]], [[matern]], [[fit-params]]."
  [kernel]
  (fn [{:keys [^double nugget ^double psill ^double range]}]
    (fn ^double [^double x]
      (if (m/zero? x) 0.0
          (->> (m// x range)
               (kernel)
               (double)
               (m/- 1.0)
               (m/* psill)
               (m/+ nugget))))))

(def ^:private semivariograms {:linear linear
                               :pentaspherical pentaspherical
                               :spherical spherical
                               :gaussian gaussian
                               :exponential exponential
                               :power power
                               :tpower power
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
  "Matheron (classical) empirical semivariogram estimator."
  ^double [_ diffs]
  (stats/mean (map m/sq diffs)))

(defn cressie-estimator
  "Cressie (robust) empirical semivariogram estimator."
  [^long n diffs]
  (-> (->> (map (fn [^double diff] (m/sqrt (m/abs diff))) diffs)
           (stats/mean))
      (m/fpow 4)
      (m// (m/+ 0.457 (m// 0.494 n) (m// 0.045 (m/* n n))))))

(defn highly-robust-estimator
  "Genton's higly robust empirical semivariogram estimator."
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

  It's the same as quantile estimator for `p` equal `0.5`."
  [_ diffs]
  (m/* 2.198109338317728 (stats/median (map m/sq diffs))))

(defn ->quantile-estimator
  "Create a quantile (Armstrong and Delfiner) empirical semivariogram estimator for given `p`."
  [^double p]
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
  "Robust M-estimator (Gunst and Hartfield) for empirical semivariogram."
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
  "Length of the diagonal of bounding box of spatial points."
  ^double [dist xs]
  (let [mn (reduce v/emn xs)
        mx (reduce v/emx xs)]
    (dist mn mx)))

(defn remove-outliers-fence
  "Removes outliers using Tukey's fences (`k`=1.5)."
  [combined ys]
  (let [[^double q1 ^double q3] (stats/quantiles ys [0.25 0.75])
        iqr (m/* 1.5 (m/- q3 q1))
        lif-thr (m/- q1 iqr)
        uif-thr (m/+ q3 iqr)]
    (filter (fn [v] (m/<= lif-thr (double (v 2)) uif-thr)) combined)))

(defn remove-outliers-mad
  "Removes outliers using median absolute deviation."
  [combined ys]
  (let [m (stats/median ys)
        s (m// (stats/median-absolute-deviation ys m) 0.6745)]
    (filter (fn [v] (m/<= (m/abs (m// (m/- (double (v 2)) m) s)) 3.0)) combined)))

(defn- maybe-remove-outliers
  [combined ys remove-outliers?]
  (cond
    (= :mad remove-outliers?) (remove-outliers-mad combined ys)
    remove-outliers? (remove-outliers-fence combined ys)
    :else combined))

(defn empirical
  "Computes the empirical (experimental) semivariogram from spatial data.

  All pairs of spatial locations are formed, filtered by the cutoff distance, and grouped
  into lag bins. For each bin, a chosen estimator computes the gamma value as half the
  average squared difference between the paired observations. The result is the primary
  input to [[fit]] and [[fit-params]] for fitting a parametric semivariogram model.

  Parameters:

  - `xs` (sequence): spatial positions of observations; each element may be a scalar or
    a vector depending on the dimensionality of the data
  - `ys` (sequence): observed values at positions `xs`; must have the same length as `xs`
  - `params` (map, optional): configuration map with the following keys:
      - `:bins` (long): number of lag bins; controls the resolution of the semivariogram;
        defaults to `15`
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
          - a custom `(fn [^long n diffs])` where `n` is the bin count and `diffs` is a
            sequence of pairwise differences
      - `:quantile` (double): quantile level used by the `:quantile` estimator;
        must be in `(0.0, 1.0)`; defaults to `0.9`
      - `:remove-outliers?`: whether to remove outliers from the raw data before binning;
        defaults to `false`; accepted values are:
          - `false` - no outlier removal
          - `true` - Tukey's fences criterion (IQR-based, k=1.5) via [[remove-outliers-fence]]
          - `:mad` - median absolute deviation criterion via [[remove-outliers-mad]]
      - `:distance` (function): a binary distance function `(fn [x1 x2])` applied to
        elements of `xs`; defaults to `fastmath.distance/euclidean`

  Returns a sequence of maps sorted by ascending lag `h`, one per non-empty bin, each containing:
  - `:n` (long): number of point pairs in the bin
  - `:h` (double): mean lag distance of pairs in the bin
  - `:gamma` (double): estimated semivariance for the bin

  See also [[fit]], [[fit-params]], [[bounding-box-diagonal]],
  [[matheron-estimator]], [[cressie-estimator]], [[highly-robust-estimator]],
  [[dowd-estimator]], [[->quantile-estimator]], [[robust-m-estimator]]."
  ([xs ys] (empirical xs ys nil))
  ([xs ys {:keys [^double cutoff ^double diagonal-den ^long bins estimator ^double quantile
                  remove-outliers? distance]
           :or {bins 15 diagonal-den 3.0 estimator :classical quantile 0.9
                remove-outliers? false distance dist/euclidean}}]
   (let [cutoff (or cutoff (/ (bounding-box-diagonal distance xs) diagonal-den))
         combined (-> (map vector (range) xs ys)
                      (maybe-remove-outliers ys remove-outliers?))
         distances (->> (for [[^long id1 x1 ^double y1] combined
                              [^long id2 x2 ^double y2] combined
                              :when (< id1 id2)
                              :let [d (double (distance x1 x2))]
                              :when (and (pos? d) (<= d cutoff))]
                          (Vec2. d (m/- y1 y2)))
                        (sort-by first))
         max-dist (.x ^Vec2 (last distances))
         splits (rest (m/slice-range 0 max-dist (inc bins)))
         estimator-fn (if (fn? estimator)
                        estimator
                        (case estimator
                          :classical matheron-estimator
                          :matheron matheron-estimator
                          :cressie cressie-estimator
                          :genton highly-robust-estimator
                          :highly-robust highly-robust-estimator
                          :dowd dowd-estimator
                          :quantile (->quantile-estimator quantile)
                          :m-robust robust-m-estimator))]
     (loop [buff []
            distances distances
            [^double c & rsplits] splits]
       (let [found (take-while (fn [^Vec2 v] (m/<= (.x v) c)) distances)]
         (if-not (seq found)
           (recur buff distances rsplits)
           (let [n (count found)
                 buff (conj buff {:n n
                                  :h (stats/mean (map first found))
                                  :gamma (m/* 0.5 (double (estimator-fn n (map second found))))})]
             (if (seq rsplits)
               (recur buff (drop n distances) rsplits)
               buff))))))))

(defn- weights-n-v1  [ns vs] (map (fn [^long n ^double v] (m// n (m/max m/EPSILON (m/abs v)))) ns vs))
(defn- weights-n-v2  [ns vs] (map (fn [^long n ^double v] (m// n (m/max m/EPSILON (m/sq v)))) ns vs))
(defn- weights-n-v3  [ns vs] (map (fn [^long n ^double v] (m// n (m/max m/EPSILON (m/cb v)))) ns vs))
(defn- est-sq  ^double [^double a ^double b] (m/sq (m/- a b)))
(defn- est-abs ^double [^double a ^double b] (m/abs (m/- a b)))
(defn- est-soft ^double [^double a ^double b] (m/* 2.0 (m/dec (m/sqrt (m/inc (est-sq a b))))))
(defn- est-huber ^double [^double a ^double b] (let [d (est-sq a b)]
                                                 (if (m/<= d 1.0) d (m/dec (m/* 2.0 (est-abs a b))))))
(defn- est-cauchy ^double [^double a ^double b] (m/log (m/inc (est-sq a b))))
(defn- est-atan ^double [^double a ^double b] (m/atan (est-sq a b)))

(defn- infer-target-args
  [defaults params]
  (->> params
       (reduce (fn [buff t]
                 (if-not (t defaults)
                   (conj buff t)
                   buff)) [])))

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
    by [[empirical]]; each map must contain `:h`, `:gamma` and `:n` keys
  - `semivariogram-model` (keyword or function): the parametric model to fit; accepted values are:
      - a keyword naming a built-in model: `:linear`, `:spherical`, `:circular`,
        `:cubic`, `:pentaspherical`, `:gaussian`, `:exponential`, `:hole`, `:power`,
        `:tpower`, `:expower`, `:cauchy`, `:rational`, `:bessel`, `:matern`,
        `:hyperspherical`, `:superspherical`, `:tplstable`
      - an RBF kernel function, which is automatically wrapped via [[rbf->variogram]]
  - `params` (map, optional): configuration map with the following keys:
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
            from the fitted model after optimisation
          - a sequence of explicit per-bin weights
          - a function `(fn ^double [^double h])` returning a scale per bin; the final
            weight is N divided by the absolute value of the function result
      - `:defaults` (map): fixed parameter values excluded from optimisation; keys are a
        subset of `:nugget`, `:psill`, `:range`, `:beta`; fitting is performed only for
        the remaining parameters. When all parameters are provided, no optimisation runs.
      - `:parameter` (double): the order or shape parameter passed to model constructors
        that support it; used as the `order` argument for `:bessel`, `:matern`,
        `:hyperspherical`, `:superspherical`, and as `H` for `:tplstable`
      - `:bounds` (map): overrides for the default parameter search bounds; keys are a
        subset of `:nugget`, `:psill`, `:range`, `:beta`, each mapping to a `[lo hi]`
        vector; default bounds are:
          - `:nugget` - `[0.0, 5 * min-gamma]`
          - `:psill` - `[0.0, 5 * max-gamma]`
          - `:range` - `[EPSILON, max-h]`
          - `:beta` - `[-30.0, 30.0]`

  Returns a map containing:
  - `:semivariogram` (function): the fitted semivariogram `(fn ^double [^double x])`
  - `:semivariogram-creator` (function): the model constructor used to build the semivariogram
  - `:params` (map): all model parameters, both fitted and fixed via `:defaults`
  - `:weights` (sequence): final per-bin weights used to compute the loss
  - `:loss` (double): the final weighted loss value after fitting

  See also [[fit]], [[empirical]], [[rbf->variogram]]."
  ([empirical-semivariogram semivariogram-model] (fit-params empirical-semivariogram semivariogram-model nil))
  ([empirical-semivariogram semivariogram-model {:keys [estimation weights defaults parameter bounds]
                                                 :or {estimation :sq weights :n}}]
   (let [semivariogram-fn (if (fn? semivariogram-model)
                            (rbf->variogram semivariogram-model)
                            (if-let [->sv (->semivariograms semivariogram-model)]
                              (if parameter (->sv parameter) (->sv))
                              (semivariograms semivariogram-model)))
         
         target-args (infer-target-args defaults (if (power-semivariograms semivariogram-model)
                                                   [:nugget :psill :range :beta]
                                                   [:nugget :psill :range]))

         hs (map :h empirical-semivariogram)
         gammas (map :gamma empirical-semivariogram)

         ns (map :n empirical-semivariogram)
         ng (map m/* ns gammas)
         
         max-h (double (last hs))  
         [^double min-gamma ^double max-gamma] (stats/extent (remove m/zero? gammas) false)

         ;; hs+ -> (h+1)/(hmax+1)
         max-h+ (m/inc max-h)
         hs+ (map (fn [^double h] (m// (m/inc h) max-h+)) hs)
         
         weights-seq (cond
                       (not weights) (repeat (count empirical-semivariogram) 1.0)
                       (sequential? weights) weights
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
                               weights))

         est-fn (case estimation
                  :abs est-abs
                  :soft est-soft
                  :huber est-huber
                  :cauchy est-cauchy
                  :atan est-atan
                  est-sq)
         
         params (if-not (seq target-args)
                  defaults
                  ;; fitting is needed
                  (let [bounds-map (merge {:nugget [0.0 (m/* 5.0 min-gamma)]
                                           :psill [0.0 (m/* 5.0 max-gamma)]
                                           :range [m/EPSILON max-h]
                                           :beta [-30.0 30.0]}
                                          bounds)
                        
                        bounds (map bounds-map target-args)

                        target (condp = weights
                                 ;; cressie
                                 :ngg (fn [& r] (let [sv (semivariogram-fn (merge (zipmap target-args r) defaults))
                                                     svh (map sv hs)
                                                     weights-seq (weights-n-v2 ns svh)]
                                                 (-> (map (fn [^double gamma- ^double gamma]
                                                            (est-fn gamma- gamma)) svh gammas)
                                                     (v/dot weights-seq))))
                                 

                                 :ngg2 (fn [& r] (let [sv (semivariogram-fn (merge (zipmap target-args r) defaults))
                                                      svh (map sv hs)
                                                      weights-seq (weights-n-v3 ng svh)]
                                                  (-> (map (fn [^double gamma- ^double gamma]
                                                             (est-fn gamma- gamma)) svh gammas)
                                                      (v/dot weights-seq))))

                                 (fn [& r] (let [sv (semivariogram-fn (merge (zipmap target-args r) defaults))]
                                            (-> (map (fn [^double h ^double gamma]
                                                       (est-fn (sv h) gamma)) hs gammas)
                                                (v/dot weights-seq)))))
                        
                        m (optim/scan-and-minimize :lbfgsb target {:bounds bounds :jitter 0.05})]
                    (->> m
                         (first)
                         (zipmap target-args)
                         (merge defaults))))
         semivariogram (semivariogram-fn params)
         final-weights (case weights-seq
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
      - an RBF kernel function, automatically wrapped via [[rbf->variogram]]
  - `parameters` (map, optional): configuration map; all keys are described in full in
    [[fit-params]]; the most commonly used keys are:
      - `:estimation` (keyword): loss function; `:sq` least squares (default), `:abs`,
        `:soft`, `:huber`, `:cauchy`, `:atan`
      - `:weights`: per-bin weighting scheme; defaults to `:n` (number of pairs); accepts
        `:n`, `:nh`, `:nh+`, `:nhh`, `:nhh+`, `:nexp`, `:nsqrt`, `:ngg`, `:ngg2`,
        a sequence of explicit weights, or a function `(fn ^double [^double h])`
      - `:defaults` (map): fixed parameter values excluded from optimisation; keys are a
        subset of `:nugget`, `:psill`, `:range`, `:beta`
      - `:parameter` (double): order for `:bessel`, `:matern`, `:hyperspherical`,
        `:superspherical`; Hurst exponent `H` for `:tplstable`
      - `:bounds` (map): overrides for default parameter search bounds; keys are a subset
        of `:nugget`, `:psill`, `:range`, `:beta`, each mapping to a `[lo hi]` vector

  Returns the fitted semivariogram function `(fn ^double [^double x])`.

  See also [[fit-params]], [[empirical]], [[rbf->variogram]]."
  ([empirical-variogram semivariogram-model] (fit empirical-variogram semivariogram-model nil))
  ([empirical-variogram semivariogram-model parameters]
   (:semivariogram (fit-params empirical-variogram semivariogram-model parameters))))
