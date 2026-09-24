(ns fastmath.stats.bins
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]
           [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.stat.descriptive DescriptiveStatistics]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- mad
  ^double [avs]
  (let [m (StatUtils/percentile avs 50.0)]
    (StatUtils/percentile (v/abs (v/shift avs (m/- m))) 50.0)))

(defn- scott-fd-helper
  "Calculate number of bins based on width of the bin `h`.

  If `h` is (near) zero, falls back to the median absolute deviation of `vvs`
  as a more robust bin-width estimate. If that is also (near) zero -- e.g.
  `vvs` is heavily tied but not fully constant, so both the primary spread
  estimate and the MAD collapse to 0 while the data's range doesn't -- bins
  of that (non-zero) range would divide by ~0, so falls back to `1` bin
  directly instead (equivalent to substituting the range itself as the bin
  width, since `ceil(range/range)` is always `1`)."
  ^long [vvs ^double h]
  (let [h (if (m/<= 0.0 h m/EPSILON) (mad vvs) h)]
    (if (m/<= 0.0 h m/EPSILON)
      1
      (let [fv (first vvs)
            ^Vec2 mm (reduce (fn [^Vec2 curr ^double v]
                               (Vec2. (m/min (.x curr) v) (m/max (.y curr) v))) (Vec2. fv fv) (rest vvs))]
        (m/max 1 (long (m/ceil (m// (m/- (.y mm) (.x mm)) h))))))))

(defn sturges
  "Estimates the number of histogram bins for a sample of size `n` using Sturges' rule.

  Formula: `ceil(log2(n)) + 1`. Assumes roughly normal data; tends to
  under-smooth (recommend too few bins) for large or non-normal samples.

  Parameters:

  - `n` (`long`): The sample size.

  Returns the estimated number of bins as a positive `long` (at least `1`,
  including when `n < 1`, since `log2` of a non-positive `n` is not finite).

  See also [[rice]], [[doane]], [[scott]], [[freedman-diaconis]]."
  ^long [^long n]
  (if (m/< n 1)
    1
    (m/max 1 (long (m/inc (m/ceil (m/log2 n)))))))

(defn rice
  "Estimates the number of histogram bins for a sample of size `n` using the Rice rule.

  Formula: `ceil(2 * cbrt(n))`. Depends only on sample size (no distributional
  assumption), generally recommending more bins than [[sturges]].

  Parameters:

  - `n` (`long`): The sample size.

  Returns the estimated number of bins as a positive `long` (at least `1`).

  See also [[sturges]], [[doane]], [[scott]], [[freedman-diaconis]], [[sqrt]],
  [[terrell-scott]]."
  ^long [^long n]
  (m/max 1 (long (m/ceil (m/* 2.0 (m/cbrt n))))))

(defn sqrt
  "Estimates the number of histogram bins for a sample of size `n` using the square-root choice.

  Formula: `ceil(sqrt(n))`. Depends only on sample size; the simplest of the
  bin-count estimators, used as the default by several spreadsheet and
  plotting tools.

  Parameters:

  - `n` (`long`): The sample size.

  Returns the estimated number of bins as a positive `long` (at least `1`,
  including when `n < 1`).

  See also [[sturges]], [[rice]], [[terrell-scott]], [[doane]], [[scott]],
  [[freedman-diaconis]]."
  ^long [^long n]
  (if (m/< n 1)
    1
    (m/max 1 (long (m/ceil (m/sqrt n))))))

(defn terrell-scott
  "Estimates the number of histogram bins for a sample of size `n` using the Terrell-Scott rule.

  Formula: `ceil(cbrt(2 * n))`. Depends only on sample size; an
  asymptotically minimal-risk rule (Terrell & Scott, 1985) related to
  [[rice]]'s formula but with the factor of `2` inside, rather than outside,
  the cube root.

  Parameters:

  - `n` (`long`): The sample size.

  Returns the estimated number of bins as a positive `long` (at least `1`,
  including when `n < 1`).

  See also [[sturges]], [[rice]], [[sqrt]], [[doane]], [[scott]],
  [[freedman-diaconis]]."
  ^long [^long n]
  (if (m/< n 1)
    1
    (m/max 1 (long (m/ceil (m/cbrt (m/* 2.0 n)))))))

(defn doane
  "Estimates the number of histogram bins for `avs` using Doane's rule.

  A refinement of [[sturges]] that accounts for the sample skewness, giving
  better results than Sturges' rule for non-normal (skewed) data. `avs` must
  be the same data `n` was counted from.

  Formula: `1 + log2(n) + log2(1 + |g1|/sigma_g1)`, where `g1` is the sample
  skewness and `sigma_g1 = sqrt(6*(n-2)/((n+1)*(n+3)))` is its estimated
  standard error.

  Parameters:

  - `avs` (`double` array): The data.
  - `n` (`long`): The sample size, `(count avs)`.

  Returns the estimated number of bins as a positive `long` (at least `1`;
  always `1` when `n < 3`, since skewness is undefined below 3 points).

  See also [[sturges]], [[rice]], [[scott]], [[freedman-diaconis]]."
  ^long [^doubles avs ^long n]
  (if (m/< n 3)
    1
    (let [stats (DescriptiveStatistics. avs)]
      (m/max 1 (long (m/ceil (m/+ (m/inc (m/log2 n))
                                  (m/log2 (m/inc (m// (m/abs (.getSkewness stats))
                                                      (m/sqrt (m// (m/* 6.0 (m/- n 2.0))
                                                                   (m/* (m/inc n) (m/+ n 3.0)))))))))))))) 

(defn scott
  "Estimates the number of histogram bins for `avs` using Scott's rule.

  Bin width: `3.5 * stddev / cbrt(n)`, then bin count is the data range
  divided by that width, rounded up. Assumes roughly normal data.

  If the estimated bin width is (near) zero (e.g. `avs` has near-zero
  variance), falls back to a more robust width estimate -- see
  [[scott-fd-helper]] for the exact fallback chain and its terminal case
  (returns `1` when even that is degenerate).

  Parameters:

  - `avs` (`double` array): The data.
  - `n` (`long`): The sample size, `(count avs)`. Must match `(count avs)`;
    `avs` must be non-empty when `n >= 1`.

  Returns the estimated number of bins as a positive `long` (at least `1`,
  including when `n < 1`, since there is then no data to derive a range from).

  See also [[sturges]], [[rice]], [[doane]], [[freedman-diaconis]]."
  ^long [^doubles avs ^long n]
  (if (m/< n 1)
    1
    (let [h (m// (m/* 3.5 (m/sqrt (StatUtils/variance avs)))
                 (m/cbrt n))]
      (scott-fd-helper avs h))))

(defn freedman-diaconis
  "Estimates the number of histogram bins for `avs` using the Freedman-Diaconis rule.

  Bin width: `2 * IQR / cbrt(n)`, then bin count is the data range divided by
  that width, rounded up. Robust to outliers (uses the interquartile range
  rather than standard deviation).

  If the estimated bin width is (near) zero (e.g. `avs` has near-zero IQR),
  falls back to a more robust width estimate -- see [[scott-fd-helper]] for
  the exact fallback chain and its terminal case (returns `1` when even that
  is degenerate).

  Parameters:

  - `avs` (`double` array): The data.
  - `n` (`long`): The sample size, `(count avs)`. Must match `(count avs)`;
    `avs` must be non-empty when `n >= 1`.

  Returns the estimated number of bins as a positive `long` (at least `1`,
  including when `n < 1`, since there is then no data to derive a range from).

  See also [[sturges]], [[rice]], [[doane]], [[scott]]."
  ^long [^doubles avs ^long n]
  (if (m/< n 1)
    1
    (let [h (m// (m/* 2.0 (m/- (StatUtils/percentile avs 75)
                               (StatUtils/percentile avs 25)))
                 (m/cbrt n))]
      (scott-fd-helper avs h))))
