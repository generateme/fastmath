(ns fastmath.stats.bins
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.grid :as grid])
  (:import [fastmath.vector Vec2]
           [fastmath.java Array]
           [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.stat.correlation Covariance]
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

(defn- iqr ^double [^doubles avs] (m/- (StatUtils/percentile avs 75)
                                       (StatUtils/percentile avs 25)))

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
    (let [h (m// (m/* 2.0 (iqr avs))
                 (m/cbrt n))]
      (scott-fd-helper avs h))))

(def ^{:private true :const true :tag 'double} normal-reference-2d
  "Normal-reference constant of a 2d histogram bin width, `(48*pi)^(1/4)` (about `3.504`)."
  (m/pow (m/* 48.0 m/PI) 0.25))

(defn- line-area
  "Square of the 1d bin width `c * spread * n^(-1/3)`: the cell area for points lying on a line."
  ^double [^double c ^double spread ^long n]
  (m/sq (m/* c spread (m/pow n (m// -1.0 3.0)))))

(defn- one-zero?
  "`true` when exactly one of `a`, `b` is zero."
  [^double a ^double b]
  (not= (m/zero? a) (m/zero? b)))

(defn scott-2d
  "Estimates a 2d histogram bin area for `xs`/`ys` using the multivariate normal-reference rule (Scott, D.W. *Multivariate Density Estimation: Theory, Practice, and Visualization*, 1992, p.82), the 2-dimensional generalization of [[scott]].

  Computed per axis: `b_x = c * stddev(xs) * n^(-1/4)`, `b_y = c * stddev(ys) * n^(-1/4)`, with `c = (48*pi)^(1/4)` (about `3.504`, the exact constant for 2 dimensions, the 1d [[scott]] uses the rounded `3.5`), then combined into the bin area `b_x * b_y`. This is the area of the best rectangular bin for bivariate normal data without correlation. Callers (e.g. `fastmath.stats/estimate-grid-size`) convert the returned area into a cell size for their grid (e.g. via `fastmath.grid/area->size`).

  Points on a line parallel to an axis (exactly one axis with zero variance) have no 2d density: the area is then the square of the 1d [[scott]] bin width of the other axis, `(3.5 * stddev * n^(-1/3))^2`.

  This is the histogram form of Scott's rule; it is unrelated to the Scott's/Silverman's rule (exponent `-1/(d+4)`) used for kernel density estimation bandwidth.

  Parameters:

  - `xs`, `ys` (`double` arrays): the two coordinate axes of the data, same length.

  Returns the estimated bin area as a `double`. Less than 2 points or two constant axes give `0.0` or `##NaN` rather than an exception; interpreting it (e.g. falling back to a usable default) is left to the caller.

  See also [[scott]], [[scott-correlated]], [[freedman-diaconis-2d]]."
  ^double [^doubles xs ^doubles ys]
  (let [n (alength xs)
        vx (StatUtils/variance xs)
        vy (StatUtils/variance ys)]
    (if (and (m/>= n 2) (one-zero? vx vy))
      (line-area 3.5 (m/sqrt (m/+ vx vy)) n)
      (let [factor (m/pow (double n) -0.25)]
        (m/* (m/* normal-reference-2d (m/sqrt vx) factor)
             (m/* normal-reference-2d (m/sqrt vy) factor))))))

;; 1 - rho^2 below this is perfect correlation up to rounding (collinear points give about 1e-16, not 0)
(def ^{:private true :const true :tag 'double} perfect-correlation-tolerance (m/* 4096.0 m/MACHINE-EPSILON2))

(defn- line-variance
  "The variance along the line when the points lie on a line (exactly one constant axis, or axes correlated up to rounding), `##NaN` otherwise (also for less than 2 points or two constant axes)."
  ^double [^doubles xs ^doubles ys]
  (if (m/< (alength xs) 2)
    ##NaN
    (let [vx (StatUtils/variance xs)
          vy (StatUtils/variance ys)]
      (cond
        (and (m/zero? vx) (m/zero? vy)) ##NaN
        (one-zero? vx vy) (m/+ vx vy)
        :else (let [cov (.covariance (Covariance.) xs ys)]
                (if (m/< (m/- 1.0 (m// (m/* cov cov) (m/* vx vy))) perfect-correlation-tolerance)
                  (m/+ vx vy)
                  ##NaN))))))

(defn- scott-line-area
  "[[scott-2d]]'s area for points on a line with the variance `variance` along it."
  ^double [^double variance ^long n]
  (line-area 3.5 (m/sqrt variance) n))

(defn scott-correlated
  "Estimates a 2d histogram bin area for `xs`/`ys` using the normal-reference rule for cells of one shape and size, taking the correlation of the axes into account.

  For bivariate normal data with standard deviations `sx`, `sy` and correlation `rho`, the cell area minimizing the asymptotic mean integrated squared error is `sqrt(12 / (n * (R_x + R_y)))`, where `R_x = 1 / (8*pi * sx^3 * sy * (1-rho^2)^(3/2))` (and `R_y` alike) are the integrated squared partial derivatives of the density, so the area is `sqrt(96*pi * sx * sy * (1-rho^2)^(3/2) / (n * (1/sx^2 + 1/sy^2)))`. Correlated data gets smaller cells than with [[scott-2d]]; for uncorrelated data with equal spreads both rules give the same area. Strongly correlated data gets very small cells: the area goes to zero as `rho` goes to `±1`.

  Grid cells (`fastmath.grid`) have one size and shape, so this rule (the best area for such cells) fits them better than [[scott-2d]] (the area of the best rectangle) when the spreads of the axes differ.

  Points on a line (exactly one constant axis, or axes correlated up to rounding, `1 - rho^2` below about `1e-12`) have no 2d density: the area is then the square of the 1d [[scott]] bin width along the line, `(3.5 * s * n^(-1/3))^2`, where `s^2 = sx^2 + sy^2` is the variance along the line.

  Parameters:

  - `xs`, `ys` (`double` arrays): the two coordinate axes of the data, same length.

  Returns the estimated bin area as a `double`, `##NaN` for less than 2 points or two constant axes (left to the caller).

  See also [[scott-2d]], [[freedman-diaconis-2d]]."
  ^double [^doubles xs ^doubles ys]
  (let [n (alength xs)
        lv (line-variance xs ys)]
    (cond
      (m/< n 2) ##NaN
      (m/valid-double? lv) (scott-line-area lv n)
      :else (let [vx (StatUtils/variance xs)
                  vy (StatUtils/variance ys)
                  cov (.covariance (Covariance.) xs ys)
                  one-rho2 (m/- 1.0 (m// (m/* cov cov) (m/* vx vy)))]
              (m/sqrt (m// (m/* 96.0 m/PI (m/sqrt (m/* vx vy)) (m/pow one-rho2 1.5))
                           (m/* (double n) (m/+ (m// vx) (m// vy)))))))))

(defn freedman-diaconis-2d
  "Estimates a 2d histogram bin area for `xs`/`ys` using the 2d Freedman-Diaconis rule: [[scott-2d]] with each standard deviation replaced by the robust estimate `IQR / 1.349`.

  Robust to outliers and heavy tails, which inflate the standard deviation. For multimodal data the interquartile range overstates the spread of each mode, so the bins are coarser than with [[scott-2d]].

  When exactly one axis has zero IQR (points on a line parallel to an axis, or an axis dominated by one value), the area is the square of the 1d [[freedman-diaconis]] bin width of the other axis, `(2 * IQR * n^(-1/3))^2`.

  Parameters:

  - `xs`, `ys` (`double` arrays): the two coordinate axes of the data, same length.

  Returns the estimated bin area as a `double`. Less than 2 points or zero IQR on both axes give `0.0` or `##NaN`, left to the caller.

  See also [[scott-2d]], [[scott-correlated]], [[freedman-diaconis]]."
  ^double [^doubles xs ^doubles ys]
  (let [n (alength xs)
        ix (iqr xs)
        iy (iqr ys)]
    (if (and (m/>= n 2) (one-zero? ix iy))
      (line-area 2.0 (m/+ ix iy) n)
      (let [factor (m/pow n -0.25)]
        (m/* (m/* normal-reference-2d (m// ix 1.349) factor)
             (m/* normal-reference-2d (m// iy 1.349) factor))))))

;; ## Data-driven estimators: cross-validation and Shimazaki-Shinomoto

(def ^{:private true :const true :tag 'long} default-candidates 200)
(def ^{:private true :const true :tag 'long} default-shifts 1)

(defn- smallest-gap
  "The smallest positive difference between values of `avs`, `##Inf` when all values are equal."
  ^double [^doubles avs]
  (let [sorted (doto (aclone avs) (java.util.Arrays/sort))
        size (alength sorted)]
    (loop [i (long 1) gap ##Inf]
      (if (m/< i size)
        (let [d (m/- (Array/aget sorted i) (Array/aget sorted (m/long-dec i)))]
          (recur (m/long-inc i) (if (m/pos? d) (m/min gap d) gap)))
        gap))))

(defn- log-spaced
  "`k` values from `lo` to `hi`, evenly spaced on the log scale (just `lo` when `k` is 1 or `lo` equals `hi`)."
  [^double lo ^double hi ^long k]
  (if (or (m/< k 2) (m/== lo hi))
    [lo]
    (let [step (m// (m/log (m// hi lo)) (m/long-dec k))]
      (mapv (fn [^long i] (m/* lo (m/exp (m/* step i)))) (range k)))))

(defn- positive-integer-option
  "The value of `option` in `opts`: `default` when missing or `nil`, throws `ex-info` unless an integer at least 1."
  ^long [opts option ^long default]
  (let [v (get opts option)]
    (cond
      (nil? v) default
      (and (int? v) (m/>= (long v) 1)) (long v)
      :else (throw (ex-info (str "Invalid " (name option) ", expected an integer, at least 1.")
                            {:option option :value v :expected "an integer, at least 1"})))))

(defn- options
  "`[shifts candidates]` from an options map, defaults for missing or `nil` values; throws `ex-info` for other values than integers at least 1."
  [opts]
  [(positive-integer-option opts :shifts default-shifts)
   (positive-integer-option opts :candidates default-candidates)])

;; ### 1d

(defn- bin-counts
  "Counts of `avs` in `nbins` bins of width `h` starting at `origin`; values beyond the last bin go to the last bin."
  ^longs [^doubles avs ^double origin ^double h ^long nbins]
  (let [counts (long-array nbins)
        last-bin (m/long-dec nbins)]
    (dotimes [i (alength avs)]
      (let [idx (long (m/floor (m// (m/- (Array/aget avs i) origin) h)))]
        (Array/inc counts (m/constrain idx 0 last-bin))))
    counts))

(defn- sum-of-squares
  ^double [^longs counts]
  (loop [i (long 0) acc 0.0]
    (if (m/< i (alength counts))
      (recur (m/long-inc i) (m/+ acc (m/sq (Array/get counts i))))
      acc)))

(deftype ^:private Sample1d [^doubles avs ^double mn ^double mx ^long shifts])

(defn- cost-1d
  "Cost of `nbins` bins over the range of the sample for `method` (`:cross-validation` or `:shimazaki-shinomoto`), averaged over the shifted origins."
  ^double [method ^Sample1d sample ^long nbins]
  (let [^doubles avs (.avs sample)
        mn (.mn sample)
        shifts (.shifts sample)
        n (double (alength avs))
        h (m// (m/- (.mx sample) mn) nbins)]
    (m// (v/sum (map (fn [^long j]
                       (let [nb (if (m/zero? j) nbins (m/long-inc nbins))
                             counts (bin-counts avs (m/- mn (m// (m/* j h) (double shifts))) h nb)
                             ss (sum-of-squares counts)]
                         (if (= method :cross-validation)
                           (m// (m/- 2.0 (m// (m/* (m/inc n) ss) (m/* n n))) h)
                           (let [kbar (m// n nb)
                                 var (m/- (m// ss nb) (m/* kbar kbar))]
                             (m// (m/- (m/* 2.0 kbar) var) (m/* h h))))))
                     (range shifts)))
         shifts)))

(defn- bins-candidates
  "Candidate bin counts, ascending: log-spaced between a quarter and twenty times the Freedman-Diaconis count, at most `n` and at most the range over twice the smallest gap."
  [^doubles avs ^long n ^double range ^long k]
  (let [fd (freedman-diaconis avs n)
        gap-limit (m/max 1.0 (m/floor (m// range (m/* 2.0 (smallest-gap avs)))))
        hi (long (m/min (m/long-mult 20 fd) n gap-limit))
        lo (m/long-min hi (m/long-max 1 (long (m/floor (m// (double fd) 4.0)))))]
    (->> (log-spaced lo hi k)
         (map m/round)
         (distinct)
         (vec))))

(defn- data-driven-bins
  ^long [method ^doubles avs opts]
  (let [n (alength avs)
        mn (StatUtils/min avs)
        mx (StatUtils/max avs)
        range (m/- mx mn)]
    (if (or (m/< n 2) (not (m/pos? range)))
      1
      (let [[^long shifts ^long k] (options opts)
            sample (Sample1d. avs mn mx shifts)]
        (first (reduce (fn [[_ ^double best-cost :as acc] ^long nbins]
                         (let [c (cost-1d method sample nbins)]
                           (if (m/< c best-cost) [nbins c] acc)))
                       [1 ##Inf] (bins-candidates avs n range k)))))))

(defn cross-validation-bins
  "Estimates the number of histogram bins for `avs` by leave-one-out cross-validation (Rudemo 1982; `stone` in numpy).

  For `N` equal bins over the range of the data (width `h`), the estimated integrated squared error of the histogram is, up to a constant, `(2 - (n+1) * sum(p_k^2)) / h`, where `p_k` is the fraction of the data in the bin `k`. The bin count with the lowest cost is returned. Unlike the normal-reference rules ([[scott]], [[freedman-diaconis]]) it adapts to the shape of the data (more bins for multimodal data), at the price of more variability and a pass over the data for every candidate.

  Candidates are about `:candidates` bin counts, evenly spaced on the log scale between a quarter and twenty times the [[freedman-diaconis]] count, and not more than `n`. For discretized data (rounded values) the cost decreases without limit once the bins are narrower than the gap between values, so the bin count is limited to the range divided by twice the smallest gap (as in the reference implementation of Shimazaki and Shinomoto). Ties go to the lower bin count.

  Parameters:

  - `avs` (`double` array): the data.
  - `opts` (map, optional):
    - `:candidates` - the number of candidate bin counts, default: `200`
    - `:shifts` - the number of bin origins the cost is averaged over (shifted by a fraction of the bin width, a smoother cost), default: `1` (no averaging)

  Returns the estimated number of bins as a positive `long`, `1` for less than 2 values or zero range. Throws `ex-info` when `:candidates` or `:shifts` is not an integer at least 1.

  See also [[shimazaki-shinomoto-bins]], [[cross-validation-2d]], [[freedman-diaconis]]."
  (^long [avs] (cross-validation-bins avs nil))
  (^long [avs opts] (data-driven-bins :cross-validation avs opts)))

(defn shimazaki-shinomoto-bins
  "Estimates the number of histogram bins for `avs` with the method of Shimazaki and Shinomoto (Neural Computation 19, 2007).

  For `N` equal bins over the range of the data (width `h`), the cost `(2 * mean - variance) / h^2` of the bin counts (the biased variance, empty bins included) estimates the integrated squared error of the histogram up to a constant. The bin count with the lowest cost is returned. Like [[cross-validation-bins]] it adapts to the shape of the data; both usually agree.

  Candidates, the limit for discretized data and the options are as in [[cross-validation-bins]].

  Parameters:

  - `avs` (`double` array): the data.
  - `opts` (map, optional): `:candidates` (default `200`) and `:shifts` (default `1`), see [[cross-validation-bins]].

  Returns the estimated number of bins as a positive `long`, `1` for less than 2 values or zero range.

  See also [[cross-validation-bins]], [[shimazaki-shinomoto-2d]]."
  (^long [avs] (shimazaki-shinomoto-bins avs nil))
  (^long [avs opts] (data-driven-bins :shimazaki-shinomoto avs opts)))

;; ### 2d

(deftype ^:private CellCounts [^double sum-of-squares ^double cells-in-range])

(def ^{:private true :const true :tag 'double} max-cell-key
  "Combined cell keys are used while the number of cells in the index range is below this (`Long/MAX_VALUE / 2`)."
  4.6e18)

(defn- run-lengths-sum-of-squares
  "Sum of squared lengths of the runs of equal values in the sorted `keys`."
  ^double [^longs keys]
  (let [n (alength keys)]
    (if (m/zero? n)
      0.0
      (loop [i (long 1) run (long 1) acc 0.0]
        (if (m/< i n)
          (if (m/== (Array/get keys i) (Array/get keys (m/long-dec i)))
            (recur (m/long-inc i) (m/long-inc run) acc)
            (recur (m/long-inc i) 1 (m/+ acc (m/* (double run) (double run)))))
          (m/+ acc (m/* (double run) (double run))))))))

(defn- cell-counts
  "Occupied-cell statistics of `xs`/`ys` on the grid `g`: the sum of squared counts and the number of cells in the index range of the occupied cells.

  Cells are counted by sorting combined cell keys; when the index range is too large for a `long` key (far outliers), by a map of cell coordinates."
  ^CellCounts [g ^doubles xs ^doubles ys]
  (let [n (alength xs)
        qs (long-array n)
        rs (long-array n)]
    (dotimes [i n]
      (let [^Vec2 c (grid/coords->cell g (Array/aget xs i) (Array/aget ys i))]
        (Array/set qs i (long (.x c)))
        (Array/set rs i (long (.y c)))))
    (let [qmin (long (reduce m/long-min Long/MAX_VALUE qs))
          qmax (long (reduce m/long-max Long/MIN_VALUE qs))
          rmin (long (reduce m/long-min Long/MAX_VALUE rs))
          rmax (long (reduce m/long-max Long/MIN_VALUE rs))
          q-span (m/inc (m/- (double qmax) (double qmin)))
          r-span (m/inc (m/- (double rmax) (double rmin)))
          cells (m/* q-span r-span)
          ss (if (m/< cells max-cell-key)
               (let [width (long r-span)
                     keys (long-array n)]
                 (dotimes [i n]
                   (Array/set keys i (m/long-add (m/long-mult (m/long-sub (Array/get qs i) qmin) width)
                                                 (m/long-sub (Array/get rs i) rmin))))
                 (java.util.Arrays/sort keys)
                 (run-lengths-sum-of-squares keys))
               (reduce (fn [^double acc ^long c] (m/+ acc (m/* (double c) (double c))))
                       0.0 (vals (frequencies (map vector qs rs)))))]
      (CellCounts. ss cells))))

(deftype ^:private Sample2d [^doubles xs ^doubles ys type ^double svx ^double svy ^long shifts])

(defn- cost-2d
  "Cost of cells of `size` for `method`, averaged over the shifted grid offsets."
  ^double [method ^Sample2d sample ^double size]
  (let [^doubles xs (.xs sample)
        ^doubles ys (.ys sample)
        type (.type sample)
        svx (.svx sample)
        svy (.svy sample)
        shifts (.shifts sample)
        n (double (alength xs))]
    (m// (double (reduce m/+ 0.0
                         (map (fn [^long j]
                                (let [shift (m// (m/* j size) (double shifts))
                                      g (grid/grid type size [(m/+ svx shift) (m/+ svy shift)])
                                      a (grid/area g)
                                      ^CellCounts cc (cell-counts g xs ys)
                                      ss (.sum-of-squares cc)]
                                  (if (= method :cross-validation)
                                    (m// (m/- 2.0 (m// (m/* (m/inc n) ss) (m/* n n))) a)
                                    (let [cells (.cells-in-range cc)
                                          kbar (m// n cells)
                                          var (m/- (m// ss cells) (m/* kbar kbar))]
                                      (m// (m/- (m/* 2.0 kbar) var) (m/* a a))))))
                              (range shifts))))
         (double shifts))))

(defn- grid-candidates
  "Candidate cell sizes, descending: log-spaced between a twentieth and four times the size of the [[freedman-diaconis-2d]] area (or the [[scott-2d]] area), not below the size of the area `4 * gap_x * gap_y`."
  [^doubles xs ^doubles ys type ^long k]
  (let [usable (fn [^double a] (when (and (m/valid-double? a) (m/pos? a)) a))
        base-area (or (usable (freedman-diaconis-2d xs ys)) (usable (scott-2d xs ys)))
        min-area (m/* 4.0 (smallest-gap xs) (smallest-gap ys))]
    (when (and base-area (m/valid-double? min-area))
      (let [base (grid/area->size type base-area)
            min-size (grid/area->size type min-area)
            sizes (filterv (fn [^double s] (m/>= s min-size))
                           (log-spaced (m// base 20.0) (m/* base 4.0) k))]
        (rseq (if (seq sizes) sizes [min-size]))))))

(defn- data-driven-area
  [method ^doubles xs ^doubles ys type opts]
  (let [[^long shifts ^long k] (options opts)
        [svx svy] (or (:sv opts) [0.0 0.0])
        n (alength xs)
        vx (StatUtils/variance xs)
        vy (StatUtils/variance ys)
        ;; a constant axis leaves no cells to compare (no gap along it), but a diagonal line or a
        ;; far outlier (correlation of 1 up to rounding) is searched as any other data
        axis-line? (and (m/>= n 2) (one-zero? vx vy))
        sizes (when (and (m/>= n 2) (not axis-line?)) (grid-candidates xs ys type k))
        sample (Sample2d. xs ys type (double svx) (double svy) shifts)]
    (cond
      axis-line? (scott-line-area (m/+ vx vy) n)
      (not (seq sizes)) ##NaN
      :else
      (let [best (first (reduce (fn [[_ ^double best-cost :as acc] ^double size]
                                  (let [c (cost-2d method sample size)]
                                    (if (m/< c best-cost) [size c] acc)))
                                [(first sizes) ##Inf] sizes))]
        (grid/area (grid/grid type best))))))

(defn cross-validation-2d
  "Estimates a 2d histogram cell area for `xs`/`ys` by leave-one-out cross-validation (Rudemo 1982), for cells of the grid `type` (`fastmath.grid`).

  For cells of area `A`, the estimated integrated squared error of the histogram is, up to a constant, `(2 - (n+1) * sum(p_k^2)) / A`, where `p_k` is the fraction of the points in the cell `k`. The cell size with the lowest cost is chosen, counting the points on the same grid as `fastmath.stats/histogram2d` (cell type `type`, offset `:sv`). Unlike the normal-reference rules ([[scott-2d]], [[scott-correlated]], [[freedman-diaconis-2d]]) it adapts to the shape of the data (smaller cells for multimodal or curved data), at the price of more variability and a pass over the data for every candidate.

  Candidates are `:candidates` cell sizes, evenly spaced on the log scale between a twentieth and four times the size of the [[freedman-diaconis-2d]] area (the [[scott-2d]] area when that one is not usable). For discretized data (rounded values) the cost decreases without limit once the cells are smaller than the gaps between values, so the cell area is at least `4 * gap_x * gap_y` (the smallest positive gaps along the axes). Ties go to the larger cell.

  Parameters:

  - `xs`, `ys` (`double` arrays): the two coordinate axes of the data, same length.
  - `type` (keyword): the cell type, one of `fastmath.grid/cell-names`.
  - `opts` (map, optional):
    - `:sv` - the grid offset (translating vector), default: `[0.0 0.0]`
    - `:candidates` - the number of candidate cell sizes, default: `200`
    - `:shifts` - the number of grid offsets the cost is averaged over (shifted along the diagonal by a fraction of the size, a smoother cost), default: `1` (no averaging)

  Points on a line parallel to an axis (one constant axis) leave no cells to compare: the area is then the square of the 1d [[scott]] bin width along the line, as in [[scott-correlated]]. Other lines, and data with a far outlier, are searched as any other data.

  Throws `ex-info` when `:candidates` or `:shifts` is not an integer at least 1.

  Returns the estimated cell area as a `double`, `##NaN` for less than 2 points or two constant axes.

  See also [[shimazaki-shinomoto-2d]], [[cross-validation-bins]]."
  (^double [xs ys type] (cross-validation-2d xs ys type nil))
  (^double [xs ys type opts] (data-driven-area :cross-validation xs ys type opts)))

(defn shimazaki-shinomoto-2d
  "Estimates a 2d histogram cell area for `xs`/`ys` with the method of Shimazaki and Shinomoto (Neural Computation 19, 2007), for cells of the grid `type` (`fastmath.grid`).

  For cells of area `A`, the cost `(2 * mean - variance) / A^2` of the cell counts (the biased variance) estimates the integrated squared error of the histogram up to a constant. The counts include the empty cells in the index range of the occupied cells (all cells of the rows and columns of the grid that contain points). Like [[cross-validation-2d]] it adapts to the shape of the data; both usually agree.

  Candidates, the limit for discretized data, points on a line parallel to an axis and the options are as in [[cross-validation-2d]].

  Parameters:

  - `xs`, `ys` (`double` arrays): the two coordinate axes of the data, same length.
  - `type` (keyword): the cell type, one of `fastmath.grid/cell-names`.
  - `opts` (map, optional): `:sv`, `:candidates` and `:shifts`, see [[cross-validation-2d]].

  Returns the estimated cell area as a `double`, `##NaN` for less than 2 points or two constant axes.

  See also [[cross-validation-2d]], [[shimazaki-shinomoto-bins]]."
  (^double [xs ys type] (shimazaki-shinomoto-2d xs ys type nil))
  (^double [xs ys type opts] (data-driven-area :shimazaki-shinomoto xs ys type opts)))
