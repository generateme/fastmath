(ns fastmath.stats
  "Statistics functions.

  * Descriptive statistics.
  * Correlation / covariance
  * Outliers
  * Confidence intervals
  * Extents
  * Effect size
  * Tests
  * Histogram
  * ACF/PACF
  * Bootstrap
  * Binary measures

  All functions are backed by Apache Commons Math or SMILE libraries. All work with Clojure sequences.

  ### Descriptive statistics

  All in one function [[stats-map]] contains:

  * `:Size` - size of the samples, `(count ...)`
  * `:Min` - [[minimum]] value
  * `:Max` - [[maximum]] value
  * `:Range` - range of values
  * `:Mean` - [[mean]]/average
  * `:Median` - [[median]], see also: [[median-3]]
  * `:Mode` - [[mode]], see also: [[modes]]
  * `:Q1` - first quartile, use: [[percentile]], [[quartile]]
  * `:Q3` - third quartile, use: [[percentile]], [[quartile]]
  * `:Total` - [[sum]] of all samples
  * `:SD` - sample standard deviation
  * `:Variance` - variance
  * `:MAD` - [[median-absolute-deviation]]
  * `:SEM` - standard error of mean
  * `:LAV` - lower adjacent value, use: [[adjacent-values]]
  * `:UAV` - upper adjacent value, use: [[adjacent-values]]
  * `:IQR` - interquartile range, `(- q3 q1)`
  * `:LOF` - lower outer fence, `(- q1 (* 3.0 iqr))`
  * `:UOF` - upper outer fence, `(+ q3 (* 3.0 iqr))`
  * `:LIF` - lower inner fence, `(- q1 (* 1.5 iqr))`
  * `:UIF` - upper inner fence, `(+ q3 (* 1.5 iqr))`
  * `:Outliers` - list of [[outliers]], samples which are outside outer fences
  * `:Kurtosis` - [[kurtosis]]
  * `:Skewness` - [[skewness]]

  Note: [[percentile]] and [[quartile]] can have 10 different interpolation strategies. See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.html)"
  {:metadoc/categories {:stat "Descriptive statistics"
                        :corr "Correlation"
                        :extent "Extents"
                        :time "Time series"
                        :effect "Effect size"
                        :test "Hypothesis test"
                        :norm "Normalize"}}
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.distance :as d]
            [fastmath.vector :as v]
            [fastmath.interpolation :as interp]
            [fastmath.optimization :as opt]
            [fastmath.kernel :as k])
  (:import [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.stat.descriptive.rank Percentile Percentile$EstimationType]
           [org.apache.commons.math3.stat.descriptive.moment Kurtosis Skewness]
           [org.apache.commons.math3.stat.correlation KendallsCorrelation SpearmansCorrelation PearsonsCorrelation]
           [org.apache.commons.math3.stat.regression SimpleRegression]
           [org.apache.commons.math3.analysis.integration RombergIntegrator]
           [org.apache.commons.math3.analysis UnivariateFunction]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(def
  ^{:metadoc/categories #{:stat}
    :docs "List of estimation strategies for [[percentile]]/[[quantile]] functions."}
  estimation-strategies-list {:legacy Percentile$EstimationType/LEGACY
                              :r1 Percentile$EstimationType/R_1
                              :r2 Percentile$EstimationType/R_2
                              :r3 Percentile$EstimationType/R_3
                              :r4 Percentile$EstimationType/R_4
                              :r5 Percentile$EstimationType/R_5
                              :r6 Percentile$EstimationType/R_6
                              :r7 Percentile$EstimationType/R_7
                              :r8 Percentile$EstimationType/R_8
                              :r9 Percentile$EstimationType/R_9})

(defn percentile
  "Calculate percentile of a `vs`.

  Percentile `p` is from range 0-100.

  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html).

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html)

  See also [[quantile]]."
  {:metadoc/categories #{:stat}}
  (^double [vs ^double p]
   (StatUtils/percentile (m/seq->double-array vs) p))
  (^double [vs ^double p estimation-strategy]
   (let [^Percentile perc (.withEstimationType ^Percentile (Percentile.) (get estimation-strategies-list estimation-strategy Percentile$EstimationType/LEGACY))]
     (.evaluate perc (m/seq->double-array vs) p))))

(defn percentiles
  "Calculate percentiles of a `vs`.

  Percentiles are sequence of values from range 0-100.

  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html).

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html)

  See also [[quantile]]."
  {:metadoc/categories #{:stat}}
  ([vs] (percentiles vs [25 50 75 100]))
  ([vs ps] (percentiles vs ps nil))
  ([vs ps estimation-strategy]
   (let [^Percentile perc (.withEstimationType ^Percentile (Percentile.) (or (estimation-strategies-list estimation-strategy) Percentile$EstimationType/LEGACY))]
     (.setData perc (m/seq->double-array vs))
     (mapv (fn [^double p] (.evaluate perc p)) ps))))

(defn quantile
  "Calculate quantile of a `vs`.

  Quantile `q` is from range 0.0-1.0.

  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html) for interpolation strategy.

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html)

  See also [[percentile]]."
  {:metadoc/categories #{:stat}}
  (^double [vs ^double q]
   (percentile vs (m/constrain (* q 100.0) 0.0 100.0)))
  (^double [vs ^double q estimation-strategy]
   (percentile vs (m/constrain (* q 100.0) 0.0 100.0) estimation-strategy)))

(defn quantiles
  "Calculate quantiles of a `vs`.

  Quantilizes is sequence with values from range 0.0-1.0.

  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html) for interpolation strategy.

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html)

  See also [[percentiles]]."
  {:metadoc/categories #{:stat}}
  ([vs] (quantiles vs [0.25 0.5 0.75 1.0]))
  ([vs qs]
   (percentiles vs (map #(m/constrain (* ^double % 100.0) 0.0 100.0) qs)))
  ([vs qs estimation-strategy]
   (percentiles vs (map #(m/constrain (* ^double % 100.0) 0.0 100.0) qs) estimation-strategy)))

(defn- wquantile-interpolator
  [vs ws method]
  (let [sorted (sort-by first (map vector vs ws))
        probabilities (map second sorted)
        data (map first sorted)
        data (conj data (first data))
        ^double wsum (reduce m/fast+ probabilities)
        weights (conj (reductions m/fast+ (map (fn [^double p] (/ p wsum)) probabilities)) 0.0)]
    (case method
      :linear (interp/linear-smile weights data)
      :average (let [interp1 (interp/step-before weights data)
                     interp2 (interp/step-after weights data)]
                 (fn [^double x] (* 0.5 (+ ^double (interp1 x)
                                          ^double (interp2 x)))))
      :step (interp/step-before weights data))))

;; based on spatstat.geom::weighted.quantile

(defn wquantile
  "Weighted quantile.

  Calculation is done using interpolation. There are three methods:
  * `:linear` - linear interpolation, default
  * `:step` - step interpolation
  * `:average` - average of ties

  Based on `spatstat.geom::weighted.quantile` from R."
  (^double [vs ws ^double q] (wquantile vs ws q :linear))
  (^double [vs ws ^double q method]
   (let [interp (wquantile-interpolator vs ws method)]
     (interp q))))

(defn wquantiles
  "Weighted quantiles.

  Calculation is done using interpolation. There are three methods:
  * `:linear` - linear interpolation, default
  * `:step` - step interpolation
  * `:average` - average of ties

  Based on `spatstat.geom::weighted.quantile` from R."
  ([vs ws] (wquantiles vs ws [0.25 0.5 0.75 1.0]))
  ([vs ws qs] (wquantiles vs ws qs :linear))
  ([vs ws qs method]
   (let [interp (wquantile-interpolator vs ws method)]
     (mapv interp qs))))

(defn wmedian
  "Weighted median.

  Calculation is done using interpolation. There are three methods:
  * `:linear` - linear interpolation, default
  * `:step` - step interpolation
  * `:average` - average of ties

  Based on `spatstat.geom::weighted.quantile` from R."
  (^double [vs ws] (wquantile vs ws 0.5))
  (^double [vs ws method] (wquantile vs ws 0.5 method)))

(defn median
  "Calculate median of `vs`. See [[median-3]]."
  {:metadoc/categories #{:stat}}
  (^double [vs estimation-strategy]
   (percentile vs 50.0 estimation-strategy))
  (^double [vs]
   (percentile vs 50.0)))

(defn median-3
  "Median of three values. See [[median]]."
  {:metadoc/categories #{:stat}}
  ^double [^double a ^double b ^double c]
  (m/max (m/min a b) (m/min (m/max a b) c)))

(defn mean
  "Calculate mean of `vs`"
  {:metadoc/categories #{:stat}}
  (^double [vs] (StatUtils/mean (m/seq->double-array vs))))

(defn geomean
  "Geometric mean for positive values only"
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (m/exp (mean (map (fn [^double v] (m/log v)) vs))))

(defn harmean
  "Harmonic mean"
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (/ (mean (map (fn [^double v] (/ v)) vs))))

(defn powmean
  "Generalized power mean"
  {:metadoc/categories #{:stat}}
  ^double [vs ^double power]
  (cond
    (zero? power) (geomean vs)
    (m/one? power) (mean vs)
    (== power m/THIRD) (m/cb (mean (map #(m/cbrt %) vs)))
    (== power 0.5) (m/sq (mean (map #(m/sqrt %) vs)))
    (== power 2.0) (m/sqrt (mean (map m/sq vs)))
    (== power 3.0) (m/cbrt (mean (map m/cb vs))) 
    :else (m/pow (mean (map (fn [^double v] (m/pow v power)) vs)) (/ power))))

(defn wmean
  "Weighted mean"
  {:metadoc/categories #{:stat}}
  (^double [vs] (mean vs))
  (^double [vs weights]
   (/ ^double (reduce m/fast+ (map m/fast* vs weights))
      ^double (reduce m/fast+ weights))))

(defn population-variance
  "Calculate population variance of `vs`.

  See [[variance]]."
  {:metadoc/categories #{:stat}}
  (^double [vs]
   (StatUtils/populationVariance (m/seq->double-array vs)))
  (^double [vs ^double u]
   (StatUtils/populationVariance (m/seq->double-array vs) u)))

(defn variance
  "Calculate variance of `vs`.

  See [[population-variance]]."
  {:metadoc/categories #{:stat}}
  (^double [vs]
   (StatUtils/variance (m/seq->double-array vs)))
  (^double [vs ^double u]
   (StatUtils/variance (m/seq->double-array vs) u)))

(defn population-stddev
  "Calculate population standard deviation of `vs`.

  See [[stddev]]."
  {:metadoc/categories #{:stat}}
  (^double [vs]
   (m/sqrt (population-variance vs)))
  (^double [vs u]
   (m/sqrt (population-variance vs u))))

(defn stddev
  "Calculate standard deviation of `vs`.

  See [[population-stddev]]."
  {:metadoc/categories #{:stat}}
  (^double [vs]
   (m/sqrt (variance vs)))
  (^double [vs u]
   (m/sqrt (variance vs u))))

(defn variation
  "Coefficient of variation CV = stddev / mean"
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [vs (m/seq->double-array vs)]
    (/ (stddev vs)
       (mean vs))))

(defn median-absolute-deviation
  "Calculate MAD"
  {:metadoc/categories #{:stat}}
  (^double [vs] (median-absolute-deviation vs nil))
  (^double [vs center]
   (let [m (double (or center (median vs)))]
     (median (map (fn [^double x] (m/abs (- x m))) vs))))
  (^double [vs center estimation-strategy]
   (let [m (double (or center (median vs estimation-strategy)))]
     (median (map (fn [^double x] (m/abs (- x m))) vs) estimation-strategy))))

(def ^{:doc "Alias for [[median-absolute-deviation]]"
     :metadoc/categories #{:stat}}
  mad median-absolute-deviation)

(defn mean-absolute-deviation
  "Calculate mean absolute deviation"
  {:metadoc/categories #{:stat}}
  (^double [vs] (mean-absolute-deviation vs nil))
  (^double [vs center]
   (let [m (double (or center (mean vs)))]
     (mean (map (fn [^double x] (m/abs (- x m))) vs)))))

(defn sem
  "Standard error of mean"
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [s (stddev vs)]
    (/ s (m/sqrt (count vs)))))

(defmacro ^:private build-extent
  [nm mid ext]
  `(defn ~nm
     ~(str " -/+ " ext " and " mid)
     {:metadoc/categories #{:extent}}
     [~'vs]
     (let [vs# (m/seq->double-array ~'vs)
           m# (~mid vs#)
           s# (~ext vs#)]
       [(- m# s#) (+ m# s#) m#])))

(build-extent stddev-extent mean stddev)
(build-extent mad-extent median median-absolute-deviation)
(build-extent sem-extent mean sem)

(defn percentile-extent
  "Return percentile range and median.

  `p` - calculates extent of `p` and `100-p` (default: `p=25`)"
  {:metadoc/categories #{:extent}}
  ([vs] (percentile-extent vs 25.0))
  ([vs ^double p] (percentile-extent vs p (- 100.0 p)))
  ([vs p1 p2] (percentile-extent vs p1 p2 :legacy))
  ([vs ^double p1 ^double p2 estimation-strategy]
   (let [avs (m/seq->double-array vs)]
     [(percentile avs p1 estimation-strategy)
      (percentile avs p2 estimation-strategy)
      (median avs)])))

(defn quantile-extent
  "Return quantile range and median.

  `q` - calculates extent of `q` and `1.0-q` (default: `q=0.25`)"
  {:metadoc/categories #{:extent}}
  ([vs] (quantile-extent vs 0.25))
  ([vs ^double q] (quantile-extent vs q (- 1.0 q)))
  ([vs q1 q2] (quantile-extent vs q1 q2 :legacy))
  ([vs ^double q1 ^double q2 estimation-strategy]
   (let [avs (m/seq->double-array vs)]
     [(quantile avs q1 estimation-strategy)
      (quantile avs q2 estimation-strategy)
      (median avs)])))

(defn pi
  "Returns PI as a map, quantile intervals based on interval size.

  Quantiles are `(1-size)/2` and `1-(1-size)/2`"
  ([vs] (pi vs 0.5))
  ([vs ^double size] (pi vs size :legacy))
  ([vs ^double size estimation-strategy]
   (let [a (* 0.5 (- 1.0 size))
         avs (m/seq->double-array vs)
         q1 (quantile avs a estimation-strategy)
         q2 (quantile avs (- 1.0 a) estimation-strategy)]
     {(m/approx (* a 100.0) 2) q1
      (m/approx (* (- 1.0 a) 100.0) 2) q2})))

(defn pi-extent
  "Returns PI extent, quantile intervals based on interval size + median.

  Quantiles are `(1-size)/2` and `1-(1-size)/2`"
  ([vs] (pi-extent vs 0.5))
  ([vs ^double size] (pi-extent vs size :legacy))
  ([vs ^double size estimation-strategy]
   (let [a (* 0.5 (- 1.0 size))]
     (quantile-extent vs a (- 1.0 a) estimation-strategy))))


(defn hpdi-extent
  "Higher Posterior Density interval + median.

  `size` parameter is the target probability content of the interval."
  ([vs] (hpdi-extent vs 0.95))
  ([vs ^double size]
   (let [avs (m/seq->double-array vs)
         nsamp (alength avs)
         gap (m/constrain (m/round (* nsamp size)) 1 (dec nsamp))
         max-idx (- nsamp gap)]
     (java.util.Arrays/sort avs)
     (loop [idx (long 0)
            min-idx (long 0)
            mn Double/MAX_VALUE]
       (if (< idx max-idx)
         (let [diff (- (aget avs (+ idx gap))
                       (aget avs idx))]
           (if (< diff mn)
             (recur (inc idx) idx diff)
             (recur (inc idx) min-idx mn)))
         [(aget avs min-idx) (aget avs (+ min-idx gap)) (median avs)])))))

(defn iqr
  "Interquartile range."
  {:metadoc/categories #{:stat}}
  (^double [vs] (iqr vs :legacy))
  (^double [vs estimation-strategy]
   (let [[^double q1 ^double q3] (percentile-extent vs 25.0 75.0 estimation-strategy)]
     (- q3 q1))))

(defn adjacent-values
  "Lower and upper adjacent values (LAV and UAV).

  Let Q1 is 25-percentile and Q3 is 75-percentile. IQR is `(- Q3 Q1)`.

  * LAV is smallest value which is greater or equal to the LIF = `(- Q1 (* 1.5 IQR))`.
  * UAV is largest value which is lower or equal to the UIF = `(+ Q3 (* 1.5 IQR))`.
  * third value is a median of samples


  Optional `estimation-strategy` argument can be set to change quantile calculations estimation type. See [[estimation-strategies]]."
  {:metadoc/categories #{:extent}}
  ([vs]
   (adjacent-values vs :legacy))
  ([vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         [q1 m q3] (percentiles avs [25.0 50.0 75.0] estimation-strategy)]
     (adjacent-values avs q1 q3 m)))
  ([vs ^double q1 ^double q3 ^double m]
   (let [avs (m/seq->double-array vs)
         iqr (* 1.5 (- q3 q1))
         lav-thr (- q1 iqr)
         uav-thr (+ q3 iqr)]
     (java.util.Arrays/sort avs)
     [(first (filter #(>= (double %) lav-thr) avs))
      (last (filter #(<= (double %) uav-thr) avs))
      m])))

(defn inner-fence-extent
  "Returns LIF, UIF and median"
  ([vs] (inner-fence-extent vs :legacy))
  ([vs estimation-strategy]
   (let [[^double q1 ^double m ^double q3] (percentiles vs [25.0 50.0 75.0] estimation-strategy)
         iqr+ (* 1.5 (- q3 q1))]
     [(- q1 iqr+) (+ q3 iqr+) m])))

(defn outer-fence-extent
  "Returns LOF, UOF and median"
  ([vs] (outer-fence-extent vs :legacy))
  ([vs estimation-strategy]
   (let [[^double q1 ^double m ^double q3] (percentiles vs [25.0 50.0 75.0] estimation-strategy)
         iqr+ (* 3.0 (- q3 q1))]
     [(- q1 iqr+) (+ q3 iqr+) m])))

(defn outliers
  "Find outliers defined as values outside inner fences.

  Let Q1 is 25-percentile and Q3 is 75-percentile. IQR is `(- Q3 Q1)`.

  * LIF (Lower Inner Fence) equals `(- Q1 (* 1.5 IQR))`.
  * UIF (Upper Inner Fence) equals `(+ Q3 (* 1.5 IQR))`.

  Returns sequence.

  Optional `estimation-strategy` argument can be set to change quantile calculations estimation type. See [[estimation-strategies]]."
  {:metadoc/categories #{:stat}}
  ([vs]
   (outliers vs :legacy))
  ([vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         q1 (percentile avs 25.0 estimation-strategy)
         q3 (percentile avs 75.0 estimation-strategy)]
     (outliers avs q1 q3)))
  ([vs ^double q1 ^double q3]
   (let [ ;;avs (m/seq->double-array vs)
         iqr (* 1.5 (- q3 q1))
         lif-thr (- q1 iqr)
         uif-thr (+ q3 iqr)]
     ;; (java.util.Arrays/sort avs)
     (filter #(let [v (double %)]
                (or (< v lif-thr)
                    (> v uif-thr))) vs))))

(declare histogram)

(defn modes
  "Find the values that appears most often in a dataset `vs`.

  Returns sequence with all most appearing values in increasing order.

  See also [[mode]]."
  ([vs method] (modes vs method {}))
  ([vs method opts]
   (let [avs (m/seq->double-array vs)]
     (case method
       :histogram (let [{:keys [bins ^double step]} (histogram avs (get opts :bins :rice))
                        ibins (vec (map-indexed #(conj %2 %1) bins))]
                    (->> ibins
                         (map (fn [[^double L ^long fm ^long id]]
                                (let [^long f1 (second (get ibins (dec id) [0 0]))
                                      ^long f2 (second (get ibins (inc id) [0 0]))]
                                  [(+ L (* step (/ (- fm f1)
                                                   (- (* 2.0 fm) f1 f2)))) (- fm)])))
                         (sort-by second)
                         (map first)))
       :pearson (let [mu (mean avs)
                      m (median avs (:estimation-strategy opts))]
                  [(- (* 3.0 m) (* 2.0 mu))])
       :kde (let [kde (k/kernel-density (get opts :kernel :gaussian) avs (:bandwidth opts))]
              (->> (map (fn [^double v] [v (- ^double (kde v))]) vs)
                   (sort-by second)
                   (map first)))
       (seq ^doubles (StatUtils/mode avs)))))
  ([vs]
   (seq ^doubles (StatUtils/mode (m/seq->double-array vs)))))

(defn mode
  "Find the value that appears most often in a dataset `vs`.

  For sample from continuous distribution, three algorithms are possible:
  * `:histogram` - calculated from [[histogram]]
  * `:kde` - calculated from KDE
  * `:pearson` - mode = mean-3(median-mean)
  * `:default` - discrete mode

  Histogram accepts optional `:bins` (see [[histogram]]). KDE method accepts `:kde` for kernel name (default `:gaussian`) and `:bandwidth` (auto). Pearson can accept `:estimation-strategy` for median.

  See also [[modes]]."
  (^double [vs method] (mode vs method {}))
  (^double [vs method opts]
   (first (modes vs method opts)))
  (^double [vs]
   (let [m (StatUtils/mode (m/seq->double-array vs))]
     (aget ^doubles m 0))))

(defn minimum
  "Minimum value from sequence."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (smile.math.MathEx/min ^doubles vs)
    (reduce clojure.core/min vs)))

(defn maximum
  "Maximum value from sequence."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (smile.math.MathEx/max ^doubles vs)
    (reduce clojure.core/max vs)))

(defn span
  "Width of the sample, maximum value minus minimum value"
  ^double [vs]
  (let [avs (m/seq->double-array vs)]
    (- (maximum avs) (minimum avs))))

(defn extent
  "Return extent (min, max, mean) values from sequence"
  {:metadoc/categories #{:extent}}
  [vs]
  (let [^double fv (first vs)]
    (conj (reduce (fn [[^double mn ^double mx] ^double v]
                    [(min mn v) (max mx v)]) [fv fv] (rest vs)) (mean vs))))

(defn sum
  "Sum of all `vs` values."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (smile.math.MathEx/sum ^doubles vs)
    (reduce (fn [^double x ^double y] (+ x y)) 0.0 vs)))

(defn moment
  "Calculate moment (central or/and absolute) of given order (default: 2).

  Additional parameters as a map:

  * `:absolute?` - calculate sum as absolute values (default: `false`)
  * `:mean?` - returns mean (proper moment) or just sum of differences (default: `true`)
  * `:center` - value of center (default: `nil` = mean)
  * `:normalize?` - apply normalization by standard deviation to the order power"
  (^double [vs] (moment vs 2.0 nil))
  (^double [vs ^double order] (moment vs order nil))
  (^double [vs ^double order {:keys [absolute? center mean? normalize?]
                              :or {mean? true}}]
   (let [in (m/seq->double-array vs)
         cin (alength in)
         out (double-array cin)
         nf (if normalize? (m/pow (variance in) (* 0.5 order)) 1.0)
         ^double center (or center (mean in))
         f (cond
             (m/one? order) m/fast-identity
             (== order 2.0) m/sq
             (== order 3.0) m/cb
             (== order 4.0) (fn ^double [^double diff] (m/sq (m/sq diff)))
             :else (fn ^double [^double diff] (m/pow diff order)))
         a (if absolute? m/abs m/fast-identity)]
     (loop [idx (int 0)]
       (when (< idx cin)
         (aset out idx ^double (f (a (- (aget in idx) center))))
         (recur (inc idx))))
     (/ (if mean? (mean out) (sum out)) nf))))

(def ^{:deprecated "Use [[moment]] function"} second-moment moment)

(defn- yule-skewness
  ^double [vs ^double u]
  (let [[^double q1 ^double q2 ^double q3] (quantiles vs [u 0.5 (- 1.0 u)])]
    (/ (+ q3 (* -2.0 q2) q1)
       (- q3 q1))))

(defn skewness
  "Calculate skewness from sequence.

  Possible types: `:G1` (default), `:g1` (`:pearson`), `:b1`, `:B1` (`:yule`), `:B3`, `:skew`, `:mode` or `:median`."
  {:metadoc/categories #{:stat}}
  (^double [vs] (skewness vs :G1))
  (^double [vs typ]
   (let [vs (m/seq->double-array vs)]
     (cond
       (sequential? typ) (cond
                           (= :mode (first typ)) (let [[_ method opts] typ]
                                                   (/ (- (mean vs) (mode vs method opts)) (stddev vs)))
                           (#{:B1 :yule} (first typ)) (yule-skewness vs (second typ)))
       (= :mode typ) (/ (- (mean vs) (mode vs)) (stddev vs))
       (= :median typ) (/ (* 3.0 (- (mean vs) (median vs))) (stddev vs))
       (#{:B1 :yule} typ) (yule-skewness vs 0.25)
       (= :B3 typ) (let [v (median vs)]
                     (/ (- (mean vs) v)
                        (moment vs 1.0 {:absolute? true :center v})))
       :else (let [^Skewness k (Skewness.)
                   n (alength vs)
                   v (.evaluate k vs)]
               (cond
                 (= :b1 typ) (* v (/ (* (- n 2.0) (dec n)) (* n n)))
                 (#{:pearson :g1} typ) (* v (/ (- n 2.0) (m/sqrt (* n (dec n)))))
                 (= :skew typ) (* v (/ (- n 2.0) (* n (m/sqrt (dec n))))) ;; artificial, to match BCa skew definition
                 :else v))))))

(defn kurtosis
  "Calculate kurtosis from sequence.

  Possible typs: `:G2` (default), `:g2` (or `:excess`), `:geary` or `:kurt`."
  {:metadoc/categories #{:stat}}
  (^double [vs] (kurtosis vs nil))
  (^double [vs typ]
   (let [vs (m/seq->double-array vs)
         n (alength vs)         ]
     (if (= typ :geary)
       (/ (mean-absolute-deviation vs)
          (population-stddev vs))
       (let [^Kurtosis k (Kurtosis.)
             v (.evaluate k vs)]
         (cond
           (#{:excess :g2} typ) (/ (- (/ (* v (- n 2) (- n 3)) (dec n)) 6.0)
                                   (inc n))
           (= :kurt typ) (+ 3.0 (/ (- (/ (* v (- n 2) (- n 3)) (dec n)) 6.0)
                                   (inc n)))
           :else v))))))

(defn ci
  "T-student based confidence interval for given data. Alpha value defaults to 0.98.

  Last value is mean."
  {:metadoc/categories #{:extent}}
  ([vs] (ci vs 0.98))
  ([vs ^double alpha]
   (let [vsa (m/seq->double-array vs)
         cnt (count vs)
         dist (r/distribution :t {:degrees-of-freedom (dec cnt)})
         ^double crit-val (r/icdf dist (- 1.0 (* 0.5 (- 1.0 alpha))))
         mean-ci (/ (* crit-val (stddev vsa)) (m/sqrt cnt))
         mn (mean vsa)]
     [(- mn mean-ci) (+ mn mean-ci) mn])))

;; https://ocw.mit.edu/courses/mathematics/18-05-introduction-to-probability-and-statistics-spring-2014/readings/MIT18_05S14_Reading24.pdf
(defn bootstrap-ci
  "Bootstrap method to calculate confidence interval.

  Alpha defaults to 0.98, samples to 1000.
  Last parameter is statistical function used to measure, default: [[mean]].

  Returns ci and statistical function value."
  {:metadoc/categories #{:extent}}
  ([vs] (bootstrap-ci vs 0.98))
  ([vs alpha] (bootstrap-ci vs alpha 1000))
  ([vs alpha samples] (bootstrap-ci vs alpha samples mean))
  ([vs ^double alpha ^long samples stat-fn]
   (let [vsa (m/seq->double-array vs)
         cnt (count vs)
         dist (r/distribution :enumerated-real {:data vsa})
         ^double m (stat-fn vsa)
         deltas (m/seq->double-array (repeatedly samples #(- ^double (mean (r/->seq dist cnt)) m)))
         q1 (quantile deltas alpha)
         q2 (quantile deltas (- 1.0 alpha))]
     [(- m q1) (- m q2) m])))

(defn bootstrap
  "Generate set of samples of given size from provided data.

  Default `samples` is 200, number of `size` defaults to sample size."
  ([vs] (bootstrap vs 200))
  ([vs samples] (bootstrap vs samples (count vs)))
  ([vs samples size]
   (let [dist (r/distribution :enumerated-real {:data vs})]
     (repeatedly samples #(r/->seq dist size)))))

(defn stats-map
  "Calculate several statistics of `vs` and return as map.

  Optional `estimation-strategy` argument can be set to change quantile calculations estimation type. See [[estimation-strategies]]."
  {:metadoc/categories #{:stat}}
  ([vs] (stats-map vs :legacy))
  ([vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         sz (alength avs)
         mn (smile.math.MathEx/min avs)
         mx (smile.math.MathEx/max avs)
         sm (smile.math.MathEx/sum avs)
         u (/ sm sz)
         mdn (median avs)
         q1 (percentile avs 25.0 estimation-strategy)
         q3 (percentile avs 75.0 estimation-strategy)
         iqr (- q3 q1)
         sd (stddev avs)
         mad (median-absolute-deviation avs)
         [lav uav] (adjacent-values avs q1 q3 mdn)]
     {:Size sz
      :Min mn
      :Max mx
      :Range (- mx mn)
      :Mean u
      :Median mdn
      :Mode (mode avs)
      :Q1 q1
      :Q3 q3
      :Total sm
      :SD sd
      :Variance (* sd sd)
      :MAD mad
      :SEM (/ sd (m/sqrt sz))
      :LAV lav
      :UAV uav
      :IQR iqr
      :LOF (- q1 (* 3.0 iqr))
      :UOF (+ q3 (* 3.0 iqr))
      :LIF (- q1 (* 1.5 iqr))
      :UIF (+ q3 (* 1.5 iqr))
      :Outliers (outliers avs q1 q3)
      :Kurtosis (kurtosis avs)
      :Skewness (skewness avs)})))

(defn standardize
  "Normalize samples to have mean = 0 and stddev = 1."
  {:metadoc/categories #{:norm}}
  [vs]
  (seq ^doubles (StatUtils/normalize (m/seq->double-array vs))))

(defn robust-standardize
  "Normalize samples to have median = 0 and MAD = 1.

  If `q` argument is used, scaling is done by quantile difference (Q_q, Q_(1-q)). Set 0.25 for IQR."
  ([vs]
   (let [avs (m/seq->double-array vs)
         mad (median-absolute-deviation avs)
         md (median avs)]
     (map (fn [^double x] (/ (- x md) mad)) vs)))
  ([vs ^double q]
   (let [[^double q1 ^double md ^double q2] (quantiles vs [q 0.5 (- 1.0 q)])
         diff (m/abs (- q2 q1))]
     (map (fn [^double x] (/ (- x md) diff)) vs))))

(defn demean
  "Subtract mean from sequence"
  {:metadoc/categories #{:norm}}
  [vs]
  (let [m (mean vs)]
    (map (fn [^double v]
           (- v m)) vs)))

(defn winsor
  "Return winsorized data. Trim is done by using quantiles, by default is set to 0.2."
  ([vs] (winsor vs 0.2))
  ([vs quantile] (winsor vs quantile :legacy))
  ([vs ^double quantile estimation-strategy]
   (let [[qlow qmid qhigh] (quantiles (remove m/nan? vs)
                                      [quantile 0.5 (- 1.0 quantile)] estimation-strategy)]
     (winsor vs qlow qhigh qmid)))
  ([vs ^double low ^double high nan]
   (let [[^double low ^double high] (if (< low high) [low high] [high low])]
     (map (fn [^double v]
            (if (m/nan? v)
              nan
              (m/constrain v low high))) vs))))

(defn trim
  "Return trimmed data. Trim is done by using quantiles, by default is set to 0.2."
  ([vs] (trim vs 0.2))
  ([vs quantile] (trim vs quantile :legacy))
  ([vs ^double quantile estimation-strategy]
   (let [[qlow qmid qhigh] (quantiles (remove m/nan? vs)
                                      [quantile 0.5 (- 1.0 quantile)] estimation-strategy)]
     (trim vs qlow qhigh qmid)))
  ([vs ^double low ^double high nan]
   (let [[^double low ^double high] (if (< low high) [low high] [high low])]
     (->> (filter (fn [^double v]
                    (or (m/nan? v)
                        (<= low v high))) vs)
          (map (fn [^double v] (if (m/nan? v) nan v)))))))

(defn rescale
  "Lineary rascale data to desired range, [0,1] by default"
  ([vs] (rescale vs 0.0 1.0))
  ([vs ^double low ^double high]
   (let [avs (m/seq->double-array vs)
         mn (minimum avs)
         mx (maximum avs)
         n (m/make-norm mn mx low high)]
     (map n vs))))

(defn covariance
  "Covariance of two sequences."
  {:metadoc/categories #{:corr}}
  (^double [[vs1 vs2]] (covariance vs1 vs2))
  (^double [vs1 vs2]
   (smile.math.MathEx/cov (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn correlation
  "Correlation of two sequences."
  {:metadoc/categories #{:corr}}
  (^double [[vs1 vs2]] (correlation vs1 vs2))
  (^double [vs1 vs2]
   (smile.math.MathEx/cor (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn spearman-correlation
  "Spearman's correlation of two sequences."
  {:metadoc/categories #{:corr}}
  (^double [[vs1 vs2]] (spearman-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation ^SpearmansCorrelation (SpearmansCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn pearson-correlation
  "Pearson's correlation of two sequences."
  {:metadoc/categories #{:corr}}
  (^double [[vs1 vs2]] (pearson-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation ^PearsonsCorrelation (PearsonsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn kendall-correlation
  "Kendall's correlation of two sequences."
  {:metadoc/categories #{:corr}}
  (^double [[vs1 vs2]] (kendall-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation ^KendallsCorrelation (KendallsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn kullback-leibler-divergence
  "Kullback-Leibler divergence of two sequences."
  {:metadoc/categories #{:corr}}
  (^double [[vs1 vs2]] (kullback-leibler-divergence vs1 vs2))
  (^double [vs1 vs2]
   (smile.math.MathEx/KullbackLeiblerDivergence (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn jensen-shannon-divergence
  "Jensen-Shannon divergence of two sequences."
  {:metadoc/categories #{:corr}}
  (^double [[vs1 vs2]] (jensen-shannon-divergence vs1 vs2))
  (^double [vs1 vs2]
   (smile.math.MathEx/JensenShannonDivergence (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn coefficient-matrix
  "Generate coefficient (correlation, covariance, any two arg function) matrix from seq of seqs. Row order.

  Default method: pearson-correlation"
  {:metadoc/categories #{:corr}}
  ([vss] (coefficient-matrix vss pearson-correlation))
  ([vss measure-fn] (coefficient-matrix vss measure-fn true))
  ([vss measure-fn symmetric?]
   (if symmetric?
     (let [avss (map-indexed (fn [id v] [id (m/seq->double-array v)]) vss)
           cache (atom {})]
       (for [[^long id1 ^doubles a] avss]
         (mapv (fn [[^long id2 ^doubles b]]
                 (let [key (if (< id1 id2) [id1 id2] [id2 id1])]
                   (if (contains? @cache key)
                     (@cache key)
                     (let [cov (measure-fn a b)]
                       (swap! cache assoc key cov)
                       cov)))) avss)))
     (let [avss (map m/seq->double-array vss)]
       (for [^doubles a avss]
         (mapv #(measure-fn a ^doubles %) avss))))))

(defn correlation-matrix
  "Generate correlation matrix from seq of seqs. Row order.

  Possible measures: `:pearson` (default), `:kendall`, `:spearman`, `:kullback-leibler` and `jensen-shannon`."
  {:metadoc/categories #{:corr}}
  ([vss] (correlation-matrix vss :pearson))
  ([vss measure]
   (let [measure (get {:pearson pearson-correlation
                       :kendall kendall-correlation
                       :spearman spearman-correlation
                       :kullback-leibler kullback-leibler-divergence
                       :jensen-shannon jensen-shannon-divergence} measure pearson-correlation)]
     (coefficient-matrix vss measure true))))

(defn covariance-matrix
  "Generate covariance matrix from seq of seqs. Row order."
  {:metadoc/categories #{:corr}}
  [vss] (coefficient-matrix vss covariance true))

(defn- maybe-number->seq [v] (if (number? v) (repeat v) v))

(defn me
  "Mean error"
  {:metadoc/categories #{:stat}}
  (^double [[vs1 vs2-or-val]] (me vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (mean (map m/fast- vs1 (maybe-number->seq vs2-or-val)))))

(defn mae
  "Mean absolute error"
  {:metadoc/categories #{:stat}}
  (^double [[vs1 vs2-or-val]] (mae vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (mean (map (comp m/abs m/fast-) vs1 (maybe-number->seq vs2-or-val)))))

(defn mape
  "Mean absolute percentage error"
  {:metadoc/categories #{:stat}}
  (^double [[vs1 vs2-or-val]] (mape vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (mean (map (fn [^double a ^double b]
                (m/abs (/ (- a b) a))) vs1 (maybe-number->seq vs2-or-val)))))

(defn rss
  "Residual sum of squares"
  {:metadoc/categories #{:stat}}
  (^double [[vs1 vs2-or-val]] (rss vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (sum (map (comp m/sq m/fast-) vs1 (maybe-number->seq vs2-or-val)))))

(defn r2
  "R2"
  (^double [[vs1 vs2-or-val]] (r2 vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (- 1.0 (/ (rss vs1 vs2-or-val)
             (moment vs1 2 {:mean? false})))))

(defn mse
  "Mean squared error"
  {:metadoc/categories #{:stat}}
  (^double [[vs1 vs2-or-val]] (mse vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (mean (map (comp m/sq m/fast-) vs1 (maybe-number->seq vs2-or-val)))))

(defn rmse
  "Root mean squared error"
  {:metadoc/categories #{:stat}}
  (^double [[vs1 vs2-or-val]] (rmse vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (m/sqrt (mse vs1 vs2-or-val))))

(defn count=
  "Count equal values in both seqs. Same as [[L0]]"
  {:metadoc/categories #{:stat}}
  (^long [[vs1 vs2-or-val]] (count= vs1 vs2-or-val))
  (^long [vs1 vs2-or-val]
   (count (filter (fn [^double v] (zero? v)) (map m/fast- vs1 (maybe-number->seq vs2-or-val))))))

(def ^{:doc "Count equal values in both seqs. Same as [[count==]]"} L0 count=)

(defn L1
  "Manhattan distance"
  (^double [[vs1 vs2-or-val]] (L1 vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (d/manhattan vs1 (take (count vs1) (maybe-number->seq vs2-or-val)))))

(defn L2sq
  "Squared euclidean distance"
  (^double [[vs1 vs2-or-val]] (L2sq vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (d/euclidean-sq vs1 (take (count vs1) (maybe-number->seq vs2-or-val)))))

(defn L2
  "Euclidean distance"
  (^double [[vs1 vs2-or-val]] (L2 vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (d/euclidean vs1 (take (count vs1) (maybe-number->seq vs2-or-val)))))

(defn LInf
  "Chebyshev distance"
  (^double [[vs1 vs2-or-val]] (LInf vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (d/chebyshev vs1 (take (count vs1) (maybe-number->seq vs2-or-val)))))

(defn psnr
  "Peak signal to noise, `max-value` is maximum possible value (default: max from `vs1` and `vs2`)"
  {:metadoc/categories #{:stat}}
  (^double [[vs1 vs2-or-val]] (psnr vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [mx1 (maximum vs1)
         ^double mx2 (if (number? vs2-or-val) vs2-or-val (maximum vs2-or-val))]
     (psnr vs1 vs2-or-val (max mx1 mx2))))
  (^double [vs1 vs2-or-val ^double max-value]
   (- (* 20.0 (m/log10 max-value))
      (* 10.0 (m/log10 (mse vs1 vs2-or-val))))))

;;

(defn- scott-fd-helper
  "Calculate number of bins based on width of the bin."
  ^double [vvs ^double h]
  (let [h (if (< h m/EPSILON) (median-absolute-deviation vvs) h)]
    (if (pos? h)
      (let [[^double mn ^double mx] (extent vvs)]
        (m/ceil (/ (- mx mn) h)))
      1.0)))

(defn estimate-bins
  "Estimate number of bins for histogram.

  Possible methods are: `:sqrt` `:sturges` `:rice` `:doane` `:scott` `:freedman-diaconis` (default).

  The number returned is not higher than number of samples."
  {:metadoc/categories #{:stat}}
  ([vs] (estimate-bins vs :freedman-diaconis))
  ([vs bins-or-estimate-method]
   (if-not (keyword? bins-or-estimate-method)
     (or bins-or-estimate-method (estimate-bins vs))
     (let [n (count vs)]
       (min n (int (condp = bins-or-estimate-method
                     :sqrt (m/sqrt n)
                     :sturges (inc (m/ceil (m/log2 n)))
                     :rice (m/ceil (* 2.0 (m/cbrt n)))
                     :doane (if (< n 3) 1 (+ (inc (m/log2 n))
                                             (m/log2 (inc (/ (m/abs (skewness vs))
                                                             (m/sqrt (/ (* 6.0 (- n 2.0))
                                                                        (* (inc n) (+ n 3.0)))))))))
                     :scott (let [vvs (m/seq->double-array vs)
                                  h (/ (* 3.5 (stddev vvs))
                                       (m/cbrt n))]
                              (scott-fd-helper vvs h))
                     (let [vvs (m/seq->double-array vs)
                           h (/ (* 2.0 (iqr vvs))
                                (m/cbrt n))]
                       (scott-fd-helper vvs h)))))))))

(defn histogram
  "Calculate histogram.

  Returns map with keys:

  * `:size` - number of bins
  * `:step` - distance between bins
  * `:bins` - list of pairs of range lower value and number of hits
  * `:min` - min value
  * `:max` - max value
  * `:samples` - number of used samples

  For estimation methods check [[estimate-bins]].

  If difference between min and max values is `0`, number of bins is set to 1."
  {:metadoc/categories #{:stat}}
  ([vs] (histogram vs :freedman-diaconis))
  ([vs bins-or-estimate-method] (histogram vs bins-or-estimate-method (extent vs)))
  ([vs bins-or-estimate-method [^double mn ^double mx]]
   (let [vs (filter (fn [^double v] (<= mn v mx)) vs)
         ^long bins (estimate-bins vs bins-or-estimate-method)
         diff (- mx mn)
         bins (if (zero? diff) 1 bins)
         step (/ diff bins)
         search-array (double-array (butlast (m/slice-range mn mx (inc bins))))
         buff (long-array bins)]

     (doseq [^double v vs]
       (let [b (java.util.Arrays/binarySearch ^doubles search-array v)
             ^int pos (if (neg? b) (m/abs (+ b 2)) b)]
         (fastmath.java.Array/inc ^longs buff pos)))

     {:size bins
      :step step
      :samples (count vs)
      :min mn
      :max mx
      :bins (map vector search-array buff)})))

;;

(defn- weighted-variance-average
  ^double [groups]
  (reduce (fn [^double sum xs]
            (+ sum (* (dec (count xs))
                      (variance xs)))) 0.0 groups))

(defn pooled-variance
  "Calculate pooled variance for samples and method.

  Methods:
  * `:unbiased` - sqrt of weighted average of variances (default)
  * `:biased` - biased version of `:unbiased`
  * `:avg` - sqrt of average of variances"
  (^double [groups] (pooled-variance groups :unbiased))
  (^double [groups method]
   (let [agroups (map m/seq->double-array groups)]
     (case method
       :biased (/ (weighted-variance-average agroups)
                  ^double (reduce m/fast+ (map alength agroups)))
       :avg (/ ^double (reduce m/fast+ (map variance agroups))
               (count groups))
       (/ (weighted-variance-average agroups)
          (- ^double (reduce m/fast+ (map alength agroups)) (count groups)))))))

(defn pooled-stddev
  "Calculate pooled standard deviation for samples and method"
  (^double [groups] (m/sqrt (pooled-variance groups)))
  (^double [groups method] (m/sqrt (pooled-variance groups method))))

(defn cohens-d
  "Cohen's d effect size for two groups"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (cohens-d group1 group2))
  (^double [group1 group2] (cohens-d group1 group2 :unbiased))
  (^double [group1 group2 method]
   (let [group1 (m/seq->double-array group1)
         group2 (m/seq->double-array group2)
         diff (- (mean group1) (mean group2))]
     (/ diff (pooled-stddev [group1 group2] method)))))

(defn- effect-size-correction
  ^double [^long df]
  (- 1.0 (/ 3.0 (dec (* 4.0 df)))))

(defn cohens-d-corrected
  "Cohen's d corrected for small group size"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (cohens-d-corrected group1 group2))
  (^double [group1 group2] (cohens-d-corrected group1 group2 :unbiased))
  (^double [group1 group2 method]
   (* (effect-size-correction (if (= method :biased)
                                (+ (count group1)
                                   (count group2))
                                (+ (count group1)
                                   (count group2) -2)))
      (cohens-d group1 group2 method))))

(defn hedges-g
  "Hedges's g effect size for two groups"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (hedges-g group1 group2))
  (^double [group1 group2]
   (cohens-d group1 group2 :unbiased)))

(defn hedges-g-corrected
  "Cohen's d corrected for small group size"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (hedges-g-corrected group1 group2))
  (^double [group1 group2]
   (cohens-d-corrected group1 group2 :unbiased)))

(defn hedges-g*
  "Less biased Hedges's g effect size for two groups, J term correction."
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (hedges-g* group1 group2))
  (^double [group1 group2]
   (let [df (+ (count group1) (count group2) -2)
         df2 (* 0.5 df)
         j (m/exp (- (m/log-gamma df2)
                     (m/log (m/sqrt df2))
                     (m/log-gamma (* 0.5 (dec df)))))]
     (* j (hedges-g group1 group2)))))

(defn glass-delta
  "Glass's delta effect size for two groups"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (glass-delta group1 group2))
  (^double [group1 group2]
   (let [group2 (m/seq->double-array group2)]
     (/ (- (mean group1) (mean group2)) (stddev group2)))))

(defn means-ratio
  "Means ratio"
  (^double [[group1 group2]] (means-ratio group1 group2))
  (^double [group1 group2] (means-ratio group1 group2 false))
  (^double [group1 group2 adjusted?]
   (let [ag1 (m/seq->double-array group1)
         ag2 (m/seq->double-array group2)
         m1 (mean ag1)
         m2 (mean ag2)]
     (if-not adjusted?
       (/ m1 m2)
       (let [v1 (variance ag1)
             v2 (variance ag2)
             n1 (alength ag1)
             n2 (alength ag2)
             J (* 0.5 (- (/ v1 (* n1 m1 m1))
                         (/ v2 (* n2 m2 m2))))]
         (m/exp (+ (- (m/log m1) (m/log m2)) J)))))))

(defn means-ratio-corrected
  "Bias correced means ratio"
  (^double [[group1 group2]] (means-ratio-corrected group1 group2))
  (^double [group1 group2]
   (means-ratio group1 group2 true)))

(defn cliffs-delta
  "Cliff's delta effect size for ordinal data."
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (cliffs-delta group1 group2))
  (^double [group1 group2]
   (/ ^double (reduce m/fast+ (for [a group1
                                    b group2]
                                (m/signum (compare a b))))
      (* (count group1) (count group2)))))

;;

(defn ameasure
  "Vargha-Delaney A measure for two populations a and b"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (ameasure group1 group2))
  (^double [group1 group2]
   (let [m (count group1)
         n (count group2)
         ^double r1 (reduce m/fast+ (take m (m/rank1 (concat group1 group2))))]
     (/ (- (+ r1 r1) (* m (inc m)))
        (* 2.0 m n)))))

(defn wmw-odds
  "Wilcoxon-Mann-Whitney odds"
  (^double [[group1 group2]] (wmw-odds group1 group2))
  (^double [group1 group2]
   (m/exp (m/logit (* 0.5 (inc (cliffs-delta group1 group2)))))))

(defn- integrate-kde
  [iterations kde ranges]
  (let [^RombergIntegrator integrator
        (RombergIntegrator. iterations RombergIntegrator/ROMBERG_MAX_ITERATIONS_COUNT)
        
        f (reify UnivariateFunction (value [_ x] (kde x)))]
    (map (fn [[^double mn ^double mx]]
           (.integrate integrator Integer/MAX_VALUE f mn mx)) ranges)))

(defn p-overlap
  "Overlapping index, kernel density approximation"
  (^double [[group1 group2]] (p-overlap group1 group2))
  (^double [group1 group2] (p-overlap group1 group2 {}))
  (^double [group1 group2 {:keys [kde bandwidth ^long min-iterations ^long steps]
                           :or {kde :gaussian min-iterations 3 steps 500}}]
   (let [[kde1 _ ^double h1 ^double mn1 ^double mx1] (k/kernel-density kde group1 bandwidth true)
         [kde2 _ ^double h2 ^double mn2 ^double mx2] (k/kernel-density kde group2 bandwidth true)
         h (* 2.0 (+ h1 h2))
         mn (- (m/min mn1 mn2) h)
         mx (+ (m/max mx1 mx2) h)
         iters (m/max 2 min-iterations)
         ranges (partition 2 1 (m/slice-range mn mx steps))
         i1 (integrate-kde iters kde1 ranges)
         i2 (integrate-kde iters kde2 ranges)]
     (reduce m/fast+ (map m/fast-min i1 i2)))))

(defn cohens-u2
  "Cohen's U2, the proportion of one of the groups that exceeds the same proportion in the other group."
  (^double [[group1 group2]] (cohens-u2 group1 group2))
  (^double [group1 group2] (cohens-u2 group1 group2 :legacy))
  (^double [group1 group2 estimation-strategy]
   (let [g1 (m/seq->double-array group1)
         g2 (m/seq->double-array group2)
         target (fn [^double p]
                  (let [p (m/constrain p 4.9E-324 0.9999999999999999)                        
                        q1 (quantiles g1 [p (- 1.0 p)] estimation-strategy)
                        q2 (quantiles g2 [(- 1.0 p) p] estimation-strategy)]
                    (-> (v/sub q1 q2)
                        (v/abs)
                        (v/mn))))]
     (ffirst (opt/minimize :brent target {:bounds [[0.5 1.0]] :initial [0.75]})))))

(defn cohens-u3
  "Cohen's U3, the proportion of the second group that is smaller than the median of the first group."
  (^double [[group1 group2]] (cohens-u3 group1 group2))
  (^double [group1 group2] (cohens-u3 group1 group2 :legacy))
  (^double [group1 group2 estimation-strategy]
   (let [m (median group1 estimation-strategy)]
     (/ (count (filter (fn [^double v] (< v m)) group2))
        (double (count group2))))))
;;

(defn pearson-r
  "Pearson `r` correlation coefficient"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (pearson-r group1 group2))
  (^double [group1 group2]
   (pearson-correlation group1 group2)))

(defn r2-determination
  "Coefficient of determination"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (r2-determination group1 group2))
  (^double [group1 group2]
   (m/sq (pearson-correlation group1 group2))))

(defn- local-linear-regression
  ^SimpleRegression [group1 group2]
  (let [lm (SimpleRegression. true)]
    (.addData lm (m/seq->double-double-array (map vector group1 group2)))
    (.regress lm)
    lm))

(defn eta-sq
  "R2, coefficient of determination"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (eta-sq group1 group2))
  (^double [group1 group2]
   (.getRSquare (local-linear-regression group1 group2))))

(defn omega-sq
  "Adjusted R2"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (omega-sq group1 group2))
  (^double [group1 group2]
   (let [lm (local-linear-regression group1 group2)
         mse (.getMeanSquareError lm)]
     (/ (- (.getRegressionSumSquares lm) mse)
        (+ (.getTotalSumSquares lm) mse)))))

(defn epsilon-sq
  "Less biased R2"
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (epsilon-sq group1 group2))
  (^double [group1 group2]
   (let [lm (local-linear-regression group1 group2)]
     (/ (- (.getRegressionSumSquares lm) (.getMeanSquareError lm))
        (.getTotalSumSquares lm)))))

(defn cohens-f2
  "Cohens f2, by default based on `eta-sq`.

  Possible `type` values are: `:eta` (default), `:omega` and `:epsilon`."
  {:metadoc/categories #{:effect}}
  (^double [[group1 group2]] (cohens-f2 group1 group2))
  (^double [group1 group2] (cohens-f2 group1 group2 :eta))
  (^double [group1 group2 type]
   (let [f (case type
             :omega omega-sq
             :epsilon epsilon-sq
             eta-sq)
         ^double v (f group1 group2)]
     (/ v (- 1.0 v)))))

(defn cohens-f
  "Cohens f, sqrt of Cohens f2.

  Possible `type` values are: `:eta` (default), `:omega` and `:epsilon`."
  (^double [[group1 group2]] (cohens-f group1 group2))
  (^double [group1 group2] (cohens-f group1 group2 :eta))
  (^double [group1 group2 type] (m/sqrt (cohens-f2 group1 group2 type))))

(defn cohens-q
  "Comparison of two correlations.

  Arity:

  * 2 - compare two correlation values
  * 3 - compare correlation of `group1` and `group2a` with correlation of `group1` and `group2b`
  * 4 - compare correlation of first two arguments with correlation of last two arguments"
  {:metadoc/categories #{:effect}}
  (^double [^double r1 ^double r2]
   (- (m/atanh r1) (m/atanh r2)))
  (^double [group1 group2a group2b]
   (cohens-q (pearson-correlation group1 group2a)
             (pearson-correlation group1 group2b)))
  (^double [group1a group2a group1b group2b]
   (cohens-q (pearson-correlation group1a group2a)
             (pearson-correlation group1b group2b))))

(declare kruskal-test)

(defn rank-eta-sq
  "Effect size for Kruskal-Wallis test"
  ^double [xs]
  (let [{:keys [^double stat ^long k ^long n]} (kruskal-test xs)]
    (m/max 0.0 (/ (inc (- stat k))
                  (- n k)))))

(defn rank-epsilon-sq
  "Effect size for Kruskal-Wallis test"
  ^double [xs]
  (let [{:keys [^double stat ^long n]} (kruskal-test xs)]
    (/ stat (/ (dec (* n n)) (inc n)))))

;;

(defn contingency-table
  "Returns frequencies map of tuples built from seqs."
  [& seqs]
  (if (= 1 (count seqs))
    (frequencies (first seqs))
    (frequencies (apply map vector seqs))))

(defn rows->contingency-table
  [xss]
  (->> (for [[row-id row] (map-indexed vector xss)
             [col-id ^long val] (map-indexed vector row)
             :when (not (zero? val))]
         [[row-id col-id] val])
       (into {})))

(defn- ct-marginals-sum
  [m]
  (->> (map (fn [[k v]]
              [k (reduce clojure.core/+ (map second v))]) m)
       (sort-by first)))

(defn contingency-table->marginals
  [ct]
  (let [rows (ct-marginals-sum (group-by ffirst ct))
        cols (ct-marginals-sum (group-by (comp second first) ct))
        n (reduce clojure.core/+ (vals ct))
        diag (filter (fn [[[a b]]] (= a b)) ct)]
    {:rows rows :cols cols :n n :diag diag}))

(defn- infer-ct
  [ct]
  (if (map? ct) ct (rows->contingency-table ct)))

(defn mcc
  "Matthews correlation coefficient also known as phi coefficient."
  ([group1 group2] (mcc (contingency-table group1 group2)))
  ([ct]
   (let [{:keys [diag rows cols ^long n]} (contingency-table->marginals (infer-ct ct))
         t (map second rows)
         p (map second cols)
         d (map second diag)
         s2 (* n n)]
     (/ (- (* n (sum d)) (v/dot t p))
        (* (m/sqrt (- s2 (v/dot p p)))
           (m/sqrt (- s2 (v/dot t t))))))))

(declare chisq-test)

(defn cramers-c
  "Cramer's C effect size for discrete data."
  {:metadoc/categories #{:effect}}
  (^double [group1 group2] (cramers-c (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (/ chi2 (+ n chi2))))))

(defn cramers-v
  "Cramer's V effect size for discrete data."
  {:metadoc/categories #{:effect}}
  (^double [group1 group2] (cramers-v (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long k ^long r ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (/ (/ chi2 n)
                (min (dec k) (dec r)))))))

(defn cramers-v-corrected
  "Corrected Cramer's V"
  {:metadoc/categories #{:effect}}
  (^double [group1 group2] (cramers-v-corrected (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long k ^long r ^long n]} (chisq-test (infer-ct contingency-table))
         k1 (dec k)
         r1 (dec r)
         n1 (double (dec n))
         phi2_ (max 0.0 (- (/ chi2 n) (/ (* k1 r1) n1)))
         k_ (- k (/ (* k1 k1) n1))
         r_ (- r (/ (* r1 r1) n1))]
     (m/sqrt (/ phi2_
                (min (dec k_) (dec r_)))))))

(defn cohens-w
  "Cohen's W effect size for discrete data."
  {:metadoc/categories #{:effect}}
  (^double [group1 group2] (cohens-w (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (/ chi2 n)))))

(defn tschuprows-t
  "Tschuprows T effect size for discrete data"
  {:metadoc/categories #{:effect}}
  (^double [group1 group2] (tschuprows-t (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long k ^long r ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (/ (/ chi2 n)
                (m/sqrt (* (dec k) (dec r))))))))

(defn cohens-kappa
  "Cohens kappa"
  (^double [group1 group2] (cohens-kappa (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^long n rows cols diag]} (contingency-table->marginals (infer-ct contingency-table))
         n2 (* n n)
         pe (/ (sum (map (fn [r c]
                           (* ^double (second r) ^double (second c))) rows cols)) n2)
         p0 (/ (sum (map second diag)) n)]
     (/ (- p0 pe) (- 1.0 pe)))))

(defn- weighted-kappa-equal-spacing
  ^double [^long R ^long id1 ^long id2]
  (- 1.0 (/ (m/abs (- id1 id2)) R)))

(defn- weighted-kappa-fleiss-cohen
  ^double [^long R ^long id1 ^long id2]
  (- 1.0 (/ (m/sq (- id1 id2)) (m/sq R))))

(defn weighted-kappa
  "Cohen's weighted kappa for indexed contingency table"
  (^double [contingency-table] (weighted-kappa contingency-table :equal-spacing))
  (^double [contingency-table weights]
   (let [ct (infer-ct contingency-table)
         R (int (reduce (fn [^long c [^long a ^long b]]
                          (max c a b)) 0 (keys ct)))
         {:keys [^long n rows cols]} (contingency-table->marginals ct)
         wfn (cond
               (= weights :equal-spacing) weighted-kappa-equal-spacing
               (= weights :fleiss-cohen) weighted-kappa-fleiss-cohen
               (fn? weights) weights
               :else (fn [_ id1 id2] (get weights [id1 id2] 0.0)))
         n2 (* n n)
         pe (/ (sum (for [[^long idr ^long r] rows
                          [^long idc ^long c] cols]
                      (* ^double (wfn R idr idc) r c))) n2)
         p0 (/ (sum (map (fn [[[^long idr ^long idc] ^long v]]
                           (* ^double (wfn R idr idc) v)) ct)) n)]
     (/ (- p0 pe) (- 1.0 pe)))))

;; binary classification statistics

(defn- binary-measures-all-calc
  [in]
  (let [{:keys [^double tp ^double fp ^double fn ^double tn] :as details} (merge {:tp 0.0
                                                                                 :fn 0.0
                                                                                 :fp 0.0
                                                                                 :tn 0.0} in)
        cp (+ tp fn)
        cn (+ fp tn)
        total (+ cp cn)
        pcp (+ tp fp)
        pcn (+ fn tn)
        ppv (/ tp pcp)
        npv (/ tn pcn)
        tpr (/ tp cp)
        fpr (/ fp cn)
        tnr (- 1.0 fpr)
        fnr (- 1.0 tpr)
        lr+ (/ tpr fpr)
        lr- (/ fnr tnr)
        ts (/ tp (+ tp fn fp))
        f-beta (clojure.core/fn [^double beta] (let [b2 (* beta beta)]
                                                 (* (inc b2) (/ (* ppv tpr)
                                                                (+ (* b2 ppv) tpr)))))
        f1-score (f-beta 1.0)
        mcc (/ (- (* tp tn) (* fp fn))
               (m/sqrt (* (+ tp fp)
                          (+ tp fn)
                          (+ tn fp)
                          (+ tn fn))))]
    (merge details {:cp cp :p cp
                    :cn cn :n cn
                    :pcp pcp :pp pcp
                    :pcn pcn :pn pcn
                    :total total
                    :tpr tpr
                    :recall tpr
                    :sensitivity tpr
                    :hit-rate tpr
                    :fnr fnr
                    :miss-rate fnr
                    :fpr fpr
                    :fall-out fpr
                    :tnr tnr
                    :specificity tnr
                    :selectivity tnr
                    :prevalence (/ cp total)
                    :accuracy (/ (+ tp tn) total)
                    :ba (/ (+ tpr tnr) 2.0)
                    :ppv ppv
                    :precision ppv
                    :fdr (- 1.0 ppv)
                    :npv npv
                    :for (- 1.0 npv)
                    :lr+ lr+
                    :lr- lr-
                    :dor (/ lr+ lr-)
                    :fm (m/sqrt (* ppv tpr))
                    :pt (/ (- (m/sqrt (* tpr fpr)) fpr)
                           (- tpr fpr))
                    :ts ts
                    :jaccard ts
                    :f-measure f1-score
                    :f1-score f1-score
                    :f-beta f-beta
                    :mcc mcc
                    :phi mcc
                    :bm (dec (+ tpr tnr))
                    :kappa (/ (* 2.0 (- (* tp tn) (* fp fn)))
                              (+ (* (+ tp fp)
                                    (+ fp tn))
                                 (* (+ tp fn)
                                    (+ fn tn))))
                    :mk (dec (+ ppv npv))})))

(defn- binary-confusion
  [t p]
  (cond
    (and t p) :tp
    (and t (not p)) :fn
    (and (not t) p) :fp
    :else :tn))

(defn- binary-process-list
  [xs true-value]
  (if-not true-value
    xs
    (let [f (if (sequential? true-value) (set true-value) true-value)]
      (map f xs))))

(defn- infer-confusion-matrix
  [confusion-matrix]
  (cond

    (and (map? confusion-matrix)
         (every? #{[:t :p] [:t :n] [:f :p] [:f :n]} (keys confusion-matrix)))
    (update-keys confusion-matrix (fn [[a b]] (keyword (str (name a) (name b)))))

    (and (sequential? confusion-matrix)
         (= 2 (count confusion-matrix))
         (every? sequential? confusion-matrix))
    (zipmap [:tp :fn :fp :tn] (flatten confusion-matrix))

    (and (sequential? confusion-matrix)
         (every? number? confusion-matrix))
    (zipmap [:tp :fn :fp :tn] confusion-matrix)
    
    :else confusion-matrix))

(defn ->confusion-matrix
  "Convert input to confusion matrix"
  ([tp fn fp tn] {:tp tp :fn fn :fp fp :tn tn})
  ([confusion-matrix] (infer-confusion-matrix confusion-matrix))
  ([actual prediction] (->confusion-matrix actual prediction nil))
  ([actual prediction encode-true]
   (let [truth (binary-process-list actual encode-true)
         prediction (binary-process-list prediction encode-true)]
     (frequencies (map binary-confusion truth prediction)))))

(defn binary-measures-all
  "Collection of binary measures.

  Arguments:
  * `confusion-matrix` - either map or sequence with `[:tp :fn :fp :tn]` values
  
  or

  * `actual` - list of ground truth values
  * `prediction` - list of predicted values
  * `true-value` - optional, true/false encoding, what is true in `truth` and `prediction`

  `true-value` can be one of:

  * `nil` - values are treating as booleans
  * any sequence - values from sequence will be treated as `true`
  * map - conversion will be done according to provided map (if there is no correspondin key, value is treated as `false`)
  * any predicate

  https://en.wikipedia.org/wiki/Precision_and_recall"
  {:metadoc/categories #{:stat}}
  ([tp fn fp tn] (binary-measures-all-calc (binary-measures-all {:tp tp :fn fn :fp fp :tn tn})))
  ([confusion-matrix] (binary-measures-all-calc (infer-confusion-matrix confusion-matrix)))
  ([actual prediction] (binary-measures-all actual prediction nil))
  ([actual prediction true-value]
   (let [truth (binary-process-list actual true-value)
         prediction (binary-process-list prediction true-value)]
     (binary-measures-all-calc (frequencies (map binary-confusion truth prediction))))))

(defn- cm-select-keys
  [cm]
  (select-keys cm [:tp :tn :fp :fn :accuracy :fdr :f-measure :fall-out :precision :recall :sensitivity :specificity :prevalence]))

(defn binary-measures
  "Subset of binary measures. See [[binary-measures-all]].

  Following keys are returned: `[:tp :tn :fp :fn :accuracy :fdr :f-measure :fall-out :precision :recall :sensitivity :specificity :prevalence]`"
  {:metadoc/categories #{:stat}}
  ([tp fn fp tn] (cm-select-keys (binary-measures-all tp fn fp tn)))
  ([confusion-matrix] (cm-select-keys (binary-measures-all confusion-matrix)))
  ([actual prediction] (binary-measures actual prediction nil))
  ([actual prediction true-value] (cm-select-keys (binary-measures-all actual prediction true-value))))

;; contingency 2x2

(defn- contingency-2x2-measures-calc
  [^long a ^long b ^long c ^long d]
  (let [fields [:a :b :c :d]
        a (or a 0) b (or b 0) c (or c 0) d (or d 0)
        table [a b c d]
        r1 (+ a b) r2 (+ c d)
        c1 (+ a c) c2 (+ b d)
        n (+ a b c d)
        dr1 (double r1) dr2 (double r2)
        dc1 (double c1) dc2 (double c2)
        dn (double n)
        expected [(/ (* r1 c1) dn) (/ (* r1 c2) dn) (/ (* r2 c1) dn) (/ (* r2 c2) dn)]
        proportions (v/div table n)
        chi2distr (r/distribution :chi-squared {:degrees-of-freedom 1})
        chi2 (let [diff (v/sub table expected)]
               (v/sum (v/ediv (v/emult diff diff) expected)))
        yates (let [diff (v/shift (v/abs (v/sub table expected)) -0.5)]
                (v/sum (v/ediv (v/emult diff diff) expected)))
        cmh (let [nmt (/ (* r1 c1) dn)]
              (/ (m/sq (- a nmt))
                 (* nmt (/ (* r2 c2) (* dn (dec dn))))))
        phi (/ (- (* a d) (* b c))
               (m/sqrt (* r1 r2 c1 c2)))
        kappa (/ (* 2.0 (- (* a d) (* b c)))
                 (+ (* c1 r2) (* r1 c2)))
        G (/ (- (+ a d) (+ b c)) dn)
        OR (/ (* a d) (double (* b c)))
        RD (- (/ a dr1) (/ c dr2))
        EER (/ a dr1) CER (/ c dr2) ARR (- CER EER)
        NNT (if (zero? ARR) ##Inf (/ ARR))
        RR (/ (/ a dr1) (/ c dr2))
        ad2 (+ (m/sq a) (m/sq d))
        hbc (* 0.5 (+ b c))
        pcc (m/sqrt (/ chi2 (+ chi2 dn)))]
    {:n n
     :table (zipmap fields table)
     :expected (zipmap fields expected)
     :marginals {:row1 (+ a b) :row2 (+ c d)
                 :col1 (+ a c) :col2 (+ b d)
                 :total n}
     :proportions {:table (zipmap fields proportions)
                   :rows (zipmap fields [(/ a dr1) (/ b dr1) (/ c dr2) (/ d dr2)])
                   :cols (zipmap fields [(/ a dc1) (/ b dc2) (/ c dc1) (/ d dc2)])
                   :marginals (merge (zipmap [:row1 :row2] (v/div [r1 r2] n))
                                     (zipmap [:col1 :col2] (v/div [c1 c2] n)))}
     :p-values {:chi2 (r/ccdf chi2distr chi2)
                :yates (r/ccdf chi2distr yates)
                :cochran-mantel-haenszel (r/ccdf chi2distr cmh)}
     :OR OR :lOR (m/log OR) :RR RR
     :risk {:RR RR :RRR (- 1.0 RR)
            :RD RD :ES r1 :CS r2
            :EER EER :CER CER :ARR ARR :NNT NNT
            :ARI (- ARR) :NNH (- NNT) :RRI (dec RR)
            :AFe (/ (dec RR) RR) :PFu (- 1.0 RR)}
     :SE (m/sqrt (v/sum (v/reciprocal table)))
     :measures {:chi2 chi2
                :yates yates
                :cochran-mantel-haenszel cmh
                :cohens-kappa kappa
                :yules-q (/ (dec OR) (inc OR))
                :holley-guilfords-g G
                :huberts-gamma (m/sq G)
                :youdens-j (/ (- (* a d) (* b c))
                              (* dr1 dr2))
                :yules-y (let [sad (m/sqrt (* a d))
                               sbc (m/sqrt (* b c))]
                           (/ (- sad sbc) (+ sad sbc)))
                :cramers-v (m/abs phi)
                :phi phi
                :scotts-pi (/ (- (* a d) (* hbc hbc))
                              (* (+ a hbc) (+ d hbc)))
                :cohens-h (* 2.0 (- (m/asin (m/sqrt (/ a dr1)))
                                    (m/asin (m/sqrt (/ c dr2))))) ;; effectsize R
                :PCC pcc
                :PCC-adjusted (* m/SQRT2 pcc)
                :TCC (m/cos (/ m/PI (inc (m/sqrt OR))))                
                :F1 (/ (+ a a)
                       (double (+ a a b c)))
                :bangdiwalas-b (/ ad2 (+ (* dr1 dc1) (* dr2 dc2)))
                :mcnemars-chi2 (/ (m/sq (- b c)) (+ b c))
                :gwets-ac1 (let [bc (+ b c)
                                 hbc2 (* bc hbc)]
                             (/ (- ad2 hbc2)
                                (+ ad2 hbc2 (* (+ a d) (+ b c)))))}}))

(defn contingency-2x2-measures-all
  ([^long a ^long b ^long c ^long d] (contingency-2x2-measures-calc a b c d))
  ([map-or-seq]
   (if (map? map-or-seq)
     (let [{:keys [a b c d] :or {a 0 b 0 c 0 d 0}} map-or-seq]
       (contingency-2x2-measures-all a b c d))
     (apply contingency-2x2-measures-calc map-or-seq)))
  ([[^long a ^long b] [^long c ^long d]] (contingency-2x2-measures-all a b c d)))

(defn contingency-2x2-measures
  [& args]
  (let [m (apply contingency-2x2-measures-all args)]
    (-> (:measures m)
        (assoc :OR (:OR m)))))

;; acf/pacf

(defn- cov-for-acf
  ^double [xs1 xs2]
  (reduce m/fast+ 0.0 (map m/fast* xs1 xs2)))

;; http://feldman.faculty.pstat.ucsb.edu/174-03/lectures/l12
(defn acf
  "Calculate acf (autocorrelation function) for given number of lags or a list of lags.

  If lags is omitted function returns maximum possible number of lags.

  See also [[acf-ci]], [[pacf]], [[pacf-ci]]"
  {:metadoc/categories #{:time}}
  ([data] (acf data (dec (count data))))
  ([data lags]
   (let [vdata (vec (demean data))
         rc (/ (double (count data)))
         lag0 (* rc (cov-for-acf vdata vdata))
         f (/ lag0)]
     (map (fn [^long lag]
            (if (zero? lag)
              1.0
              (let [v2 (subvec vdata lag)
                    v1 (subvec vdata 0 (count v2))]
                (* f rc (cov-for-acf v1 v2))))) (if (number? lags)
                                                  (range (inc (int lags)))
                                                  (seq lags))))))

;; http://feldman.faculty.pstat.ucsb.edu/174-03/lectures/l13
(defn pacf
  "Caluclate pacf (partial autocorrelation function) for given number of lags.

  If lags is omitted function returns maximum possible number of lags.

  `pacf` returns also lag `0` (which is `0.0`).

  See also [[acf]], [[acf-ci]], [[pacf-ci]]"
  {:metadoc/categories #{:time}}
  ([data] (pacf data (dec (count data))))
  ([data ^long lags]
   (let [acf (vec (acf data lags))
         phis (reductions (fn [curr ^long id]
                            (let [phi (/ (- ^double (acf id)
                                            ^double (reduce m/fast+
                                                            (map-indexed (fn [^long idx ^double c]
                                                                           (* c ^double (acf (dec (- id idx))))) curr)))
                                         (- 1.0
                                            ^double (reduce m/fast+
                                                            (map-indexed (fn [^long id ^double c]
                                                                           (* c ^double (acf (inc id)))) curr))))]

                              (conj (mapv (fn [^double p1 ^double p2]
                                            (- p1 (* phi p2))) curr (reverse curr)) phi))) [(acf 1)] (range 2 (inc lags)))]
     (conj (map last phis) 0.0))))

(defn- p-acf-ci-value
  ^double [data ^double alpha]
  (* (/ (m/sqrt (count data)))
     ^double (r/icdf r/default-normal (* 0.5 (inc (- 1.0 alpha))))))

(defn pacf-ci
  "[[pacf]] with added confidence interval data."
  {:metadoc/categories #{:time}}
  ([data] (pacf-ci data (dec (count data))))
  ([data lags] (pacf-ci data lags 0.05))
  ([data ^long lags ^double alpha]
   (let [pacf-data (pacf data lags)
         ci (p-acf-ci-value data alpha)]
     {:ci ci
      :pacf pacf-data})))

(defn acf-ci
  "[[acf]] with added confidence interval data.

  `:cis` contains list of calculated ci for every lag."
  {:metadoc/categories #{:time}}
  ([data] (acf-ci data (dec (count data))))
  ([data lags] (acf-ci data lags 0.05))
  ([data ^long lags ^double alpha]
   (let [acf-data (acf data lags)
         ci (p-acf-ci-value data alpha)]
     {:ci ci
      :acf acf-data
      :cis (map (fn [^double r]
                  (* ci (m/sqrt (dec (+ r r))))) (reductions (fn [^double acc ^double s]
                                                               (+ acc (* s s))) acf-data))})))

;;

(defn- estimate-acceleration
  "Estimates acceleration for BCA bootstrap confidence interval computation"
  ^double [avs]
  (/ (skewness avs :skew) -6.0))

(defn- cdf-accelerated-quantile
  ^double [^double z0 ^double z ^double a]
  (let [num (+ z0 z)
        denom (- 1.0 (* a num))]
    (->> (+ z0 (/ num denom))
         (r/cdf r/default-normal))))

(defn- empirical-cdf
  ^double [vs ^double value]
  (/ (double (count (filter (fn [^double v] (< v value)) vs))) (count vs)))

(defn- percentile-bca-common
  [avs p1 p2 m accel estimation-strategy]
  (let [^double z0 (r/icdf r/default-normal (empirical-cdf avs m))
        ^double z1 (r/icdf r/default-normal (/ ^double p1 100.0))
        ^double z2 (r/icdf r/default-normal (/ ^double p2 100.0))
        q1 (cdf-accelerated-quantile z0 z1 accel)
        q2 (cdf-accelerated-quantile z0 z2 accel)]
    [(quantile avs q1 estimation-strategy)
     (quantile avs q2 estimation-strategy)
     m]))

(defn percentile-bca-extent
  "Return bias corrected percentile range and mean for bootstrap samples. Also accounts for variance
   variations throught the accelaration parameter.
  See https://projecteuclid.org/euclid.ss/1032280214

  `p` - calculates extent of bias corrected `p` and `100-p` (default: `p=2.5`)

  Set `estimation-strategy` to `:r7` to get the same result as in R `coxed::bca`."
  {:metadoc/categories #{:extent}}
  ([vs] (percentile-bca-extent vs 2.5))
  ([vs ^double p] (percentile-bca-extent vs p (- 100.0 p)))
  ([vs p1 p2] (percentile-bca-extent vs p1 p2 :legacy))
  ([vs p1 p2 estimation-strategy]
   (let [avs (m/seq->double-array vs)
         accel (estimate-acceleration avs)]
     (percentile-bca-common avs p1 p2 (mean avs) accel estimation-strategy)))
  ([vs p1 p2 accel estimation-strategy]
   (let [avs (m/seq->double-array vs)]
     (percentile-bca-common avs p1 p2 (mean avs) accel estimation-strategy))))

(defn percentile-bc-extent
  "Return bias corrected percentile range and mean for bootstrap samples.
  See https://projecteuclid.org/euclid.ss/1032280214

  `p` - calculates extent of bias corrected `p` and `100-p` (default: `p=2.5`)

  Set `estimation-strategy` to `:r7` to get the same result as in R `coxed::bca`."
  {:metadoc/categories #{:extent}}
  ([vs] (percentile-bc-extent vs 2.5))
  ([vs ^double p] (percentile-bc-extent vs p (- 100.0 p)))
  ([vs p1 p2] (percentile-bc-extent vs p1 p2 :legacy))
  ([vs p1 p2 estimation-strategy]
   (percentile-bca-extent vs p1 p2 0.0 estimation-strategy)))

;;

(def binomial-ci-methods [:asymptotic :agresti-coull :clopper-pearson :wilson :prop.test
                        :cloglog :logit :probit :arcsine])

(defn binomial-ci
  "Return confidence interval for a binomial distribution.

  Possible methods are:
  * `:asymptotic` (normal aproximation, based on central limit theorem), default
  * `:agresti-coull`
  * `:clopper-pearson`
  * `:wilson`
  * `:prop.test` - one sample proportion test
  * `:cloglog`
  * `:logit`
  * `:probit`
  * `:arcsine`
  * `:all` - apply all methods and return a map of triplets

  Default confidence level: 0.95
  
  Returns a triple [lower ci, upper ci, mean]"
  {:metadoc/categories #{:extent}}
  ([^long number-of-successes ^long number-of-trials]
   (binomial-ci number-of-successes number-of-trials :asymptotic))
  ([^long number-of-successes ^long number-of-trials method]
   (binomial-ci number-of-successes number-of-trials method 0.95))
  ([^long number-of-successes ^long number-of-trials method ^double confidence-level]
   (let [p (/ (double number-of-successes) number-of-trials)
         alpha (- 1.0 confidence-level)
         alpha2 (* 0.5 alpha)
         ^double z (r/icdf r/default-normal (- 1.0 alpha2))
         z2 (* z z)
         x0? (zero? number-of-successes)
         xn? (== number-of-trials number-of-successes)]
     (case method
       :all (into {} (map #(vector %1 (binomial-ci number-of-successes number-of-trials % confidence-level))
                          binomial-ci-methods))
       :cloglog (let [logp (m/log p)
                      mu (m/log (- logp))
                      sd (* z (m/sqrt (-> (- 1.0 p) (/ number-of-trials) (/ p) (/ (* logp logp)))))
                      lcl (cond
                            x0? 0.0
                            xn? (m/pow alpha2 (/ 1.0 number-of-trials))
                            :else (m/exp (- (m/exp (+ mu sd)))))
                      ucl (cond
                            x0? (- 1.0 (m/pow alpha2 (/ 1.0 number-of-trials)))
                            xn? 1.0
                            :else (m/exp (- (m/exp (- mu sd)))))]
                  [lcl ucl p])
       :logit (let [logitp (- (m/log p) (m/log1p (- p)))
                    sd (* z (m/sqrt (-> (/ 1.0 number-of-trials) (/ p) (/ (- 1.0 p)))))
                    lcl (cond
                          x0? 0.0
                          xn? (m/pow alpha2 (/ 1.0 number-of-trials))
                          :else (let [lcl (m/exp (- logitp sd))]
                                  (/ lcl (inc lcl))))
                    ucl (cond
                          x0? (- 1.0 (m/pow alpha2 (/ 1.0 number-of-trials)))
                          xn? 1.0
                          :else (let [ucl (m/exp (+ logitp sd))]
                                  (/ ucl (inc ucl))))]
                [lcl ucl p])
       :probit (let [^double probitp (r/icdf r/default-normal p)
                     sd (* z (m/sqrt (-> (* p (- 1.0 p))
                                         (/ number-of-trials)
                                         (/ (m/sq ^double (r/pdf r/default-normal probitp))))))
                     ^double lcl (cond
                                   x0? 0.0
                                   xn? (m/pow alpha2 (/ 1.0 number-of-trials))
                                   :else (r/cdf r/default-normal (- probitp sd)))
                     ^double ucl (cond
                                   x0? (- 1.0 (m/pow alpha2 (/ 1.0 number-of-trials)))
                                   xn? 1.0
                                   :else (r/cdf r/default-normal (+ probitp sd)))]
                 [lcl ucl p])
       :prop.test (let [yatesn (/ (min 0.5 (m/abs (- number-of-successes (* number-of-trials 0.5))))
                                  number-of-trials)
                        nn (* 2.0 number-of-trials)
                        z22n (/ z2 nn)
                        z22n+ (inc (* 2.0 z22n))
                        z22n2n (/ z22n nn)
                        pc (- p yatesn)
                        pl (if-not (pos? pc) 0.0
                                   (/ (- (+ pc z22n) (* z (m/sqrt (+ (* pc (/ (- 1.0 pc) number-of-trials))
                                                                     z22n2n))))
                                      z22n+))
                        pc (+ p yatesn)
                        pu (if (>= pc 1.0) 1.0
                               (/ (+ pc z22n (* z (m/sqrt (+ (* pc (/ (- 1.0 pc) number-of-trials))
                                                             z22n2n))))
                                  z22n+))]
                    [pl pu p])
       :wilson (let [z2n (/ z2 number-of-trials)
                     p1 (+ p (* 0.5 z2n))
                     p2 (* z (m/sqrt (/ (+ (* p (- 1.0 p))
                                           (* 0.25 z2n)) number-of-trials)))
                     p3 (inc z2n)]
                 [(/ (- p1 p2) p3) (/ (+ p1 p2) p3) p])
       :clopper-pearson (let [diff (- number-of-trials
                                      number-of-successes)
                              lclbeta (if x0? 1.0
                                          ^double (r/icdf
                                                   (r/distribution :beta
                                                                   {:alpha (inc diff)
                                                                    :beta number-of-successes})
                                                   (- 1.0 alpha2)))
                              uclbeta (if xn? 0.0
                                          ^double (r/icdf
                                                   (r/distribution :beta
                                                                   {:alpha diff
                                                                    :beta (inc number-of-successes)})
                                                   alpha2))]
                          [(- 1.0 lclbeta) (- 1.0 uclbeta) p])
       :agresti-coull (let [x (+ number-of-successes (* 0.5 z2))
                            n (+ number-of-trials z2)
                            p' (/ x n)
                            zse (* z (m/sqrt (* p' (/ (- 1.0 p') n))))]
                        [(- p' zse) (+ p' zse) p])
       :arcsine (let [ap (m/asin (m/sqrt p))
                      zn (/ z (* 2.0 (m/sqrt number-of-trials)))]
                  [(m/sq (m/sin (max 0.0 (- ap zn))))
                   (m/sq (m/sin (min m/HALF_PI (+ ap zn)))) p])
       (let [zse (* z (m/sqrt (* p (/ (- 1.0 p) number-of-trials))))]
         [(- p zse) (+ p zse) p])))))

;; tests

;; t-test, reimplementation of R version

(defmacro ^:private sides-case
  [sides both right left]
  `(case ~sides
     (:two-sided :both) ~both
     (:one-sided-greater :right) ~right
     ~left))

(defn p-value
  "Calculate p-value for given distribution (default: N(0,1)), `stat`  and sides (one of `:two-sided`, `:one-sided-greater` or `:one-sided-less`/`:one-sided`)."
  ([^double stat] (p-value r/default-normal stat))
  ([distribution ^double stat] (p-value distribution stat :two-sided))
  ([distribution ^double stat sides]
   (let [stat2 (if (r/continuous? distribution) stat (dec stat))]
     (sides-case sides
                 (min 1.0 (* 2.0 (min (r/cdf distribution stat)
                                      (r/ccdf distribution stat2))))
                 (r/ccdf distribution stat2)
                 (r/cdf distribution stat)))))

(defn binomial-test
  "Binomial test

  * `alpha` - significance level (default: `0.05`)
  * `sides` - one of: `:two-sided` (default), `:one-sided-less` (short: `:one-sided`) or `:one-sided-greater`
  * `ci-method` - see [[binomial-ci-methods]]
  * `p` - tested probability"
  ([xs] (binomial-test xs {}))
  ([xs maybe-params]
   (if (map? maybe-params)
     (let [{:keys [true-false-conv] :as params} maybe-params
           xxs (binary-process-list xs true-false-conv)
           nos (count (filter identity xxs))
           not (count xxs)]
       (binomial-test nos not params))
     (binomial-test xs maybe-params {})))
  ([^long number-of-successes ^long number-of-trials {:keys [^double alpha ^double p ci-method sides]
                                                      :or {alpha 0.05 p 0.5
                                                           ci-method :asymptotic sides :two-sided}}]
   (let [distr (r/distribution :binomial {:trials number-of-trials :p p})]
     {:p-value (p-value distr (double number-of-successes) sides)
      :p p
      :successes number-of-successes
      :trials number-of-trials
      :alpha alpha
      :level (- 1.0 alpha)
      :test-type sides
      :stat number-of-successes
      :estimate (/ (double number-of-successes) number-of-trials)
      :confidence-interval (let [ci (partial binomial-ci number-of-successes number-of-trials ci-method)]
                             (sides-case sides
                                         (vec (butlast (ci (- 1.0 alpha))))
                                         [(first (ci (- 1.0 (* alpha 2.0)))) 1.0]
                                         [0.0 (second (ci (- 1.0 (* alpha 2.0))))]))})))

;; t/z

(defn- test-update-ci
  [^double mu ^double stderr [^double l ^double r]]
  [(+ mu (* l stderr)) (+ mu (* r stderr))])

(defn- test-pvalue-ci
  [d sides ^double stat ^double alpha]
  {:confidence-interval (sides-case sides
                                    (let [^double cint (r/icdf d (- 1.0 (* 0.5 alpha)))]
                                      [(- stat cint) (+ stat cint)])
                                    [(- stat ^double (r/icdf d (- 1.0 alpha))) ##Inf]
                                    [##-Inf (+ stat ^double (r/icdf d (- 1.0 alpha)))])
   :p-value (p-value d stat sides)})

(defn- test-one-sample
  [xs {:keys [^double alpha sides ^double mu]
       :or {alpha 0.05 sides :two-sided mu 0.0}}]
  (let [axs (m/seq->double-array xs)
        n (alength axs)
        m (mean axs)
        v (variance axs)
        stderr (m/sqrt (/ v n))]
    (assert (> stderr (* 10.0 m/MACHINE-EPSILON (m/abs m))) "Constant data, can't perform test.")
    {:n n
     :estimate m
     :mu mu
     :stat (/ (- m mu) stderr)
     :test-type sides
     :stderr stderr
     :alpha alpha
     :level (- 1.0 alpha)}))

(defn t-test-one-sample
  "One sample Student's t-test

  * `alpha` - significance level (default: `0.05`)
  * `sides` - one of: `:two-sided`, `:one-sided-less` (short: `:one-sided`) or `:one-sided-greater`
  * `mu` - mean (default: `0.0`)"
  {:metadoc/categories #{:test}}
  ([xs] (t-test-one-sample xs {}))
  ([xs m]
   (let [{:keys [^long n ^double stat test-type ^double alpha ^double mu ^double stderr]
          :as res} (test-one-sample xs m)
         df (dec n)
         pvals (-> (test-pvalue-ci (r/distribution :t {:degrees-of-freedom df}) test-type stat alpha)
                   (update :confidence-interval (partial test-update-ci mu stderr)))]
     (assoc (merge pvals res) :df df :t stat))))

(def ^{:deprecated "Use [[t-test-one-sample]]"} ttest-one-sample t-test-one-sample)

(defn z-test-one-sample
  "One sample z-test

  * `alpha` - significance level (default: `0.05`)
  * `sides` - one of: `:two-sided`, `:one-sided-less` (short: `:one-sided`) or `:one-sided-greater`
  * `mu` - mean (default: `0.0`)"
  {:metadoc/categories #{:test}}
  ([xs] (z-test-one-sample xs {}))
  ([xs m]
   (let [{:keys [^double stat test-type ^double alpha ^double mu ^double stderr]
          :as res} (test-one-sample xs m)
         pvals (-> (test-pvalue-ci r/default-normal test-type stat alpha)
                   (update :confidence-interval (partial test-update-ci mu stderr)))]
     (assoc (merge pvals res) :z stat))))

(defn- test-equal-variances
  [^double nx ^double ny ^double vx ^double vy]
  (let [df (- (+ nx ny) 2.0)
        v (/ (+ (* vx (dec nx))
                (* vy (dec ny))) df)]
    [df (m/sqrt (* v (+ (/ 1.0 nx)
                        (/ 1.0 ny))))]))

(defn- test-not-equal-variances
  [^double nx ^double ny ^double vx ^double vy]
  (let [stderrx (m/sqrt (/ vx nx))
        stderry (m/sqrt (/ vy ny))
        stderr (m/hypot-sqrt stderrx stderry)
        df (/ (m/sq (m/sq stderr))
              (+ (/ (m/sq (m/sq stderrx)) (dec nx))
                 (/ (m/sq (m/sq stderry)) (dec ny))))]
    [df stderr]))

(defn- test-two-samples-not-paired
  [xs ys {:keys [^double alpha sides ^double mu equal-variances?]
          :or {alpha 0.05 sides :two-sided mu 0.0 equal-variances? false}}]
  (let [axs (m/seq->double-array xs)
        ays (m/seq->double-array ys)
        nx (alength axs)
        ny (alength ays)
        mx (mean axs)
        my (mean ays)
        vx (variance axs)
        vy (variance ays)
        [df ^double stderr] (if equal-variances?
                              (test-equal-variances nx ny vx vy)
                              (test-not-equal-variances nx ny vx vy))]
    {:n [nx ny] :nx nx :ny ny
     :estimated-mu [mx my]
     :mu mu
     :estimate (- mx my mu)
     :stat (/ (- mx my mu) stderr)
     :sides sides
     :test-type sides
     :stderr stderr
     :alpha alpha
     :level (- 1.0 alpha)
     :df df
     :paired? false
     :equal-variances? equal-variances?}))

(defn t-test-two-samples
  "Two samples Student's t-test

  * `alpha` - significance level (default: `0.05`)
  * `sides` - one of: `:two-sided` (default), `:one-sided-less` (short: `:one-sided`) or `:one-sided-greater`
  * `mu` - mean (default: `0.0`)
  * `paired?` - unpaired or paired test, boolean (default: `false`)
  * `equal-variances?` - unequal or equal variances, boolean (default: `false`)"
  {:metadoc/categories #{:test}}
  ([xs ys] (t-test-two-samples xs ys {}))
  ([xs ys {:keys [paired? equal-variances?]
           :or {paired? false equal-variances? false}
           :as params}]
   (let [nx (count xs)
         ny (count ys)]
     (assert (or (and equal-variances? (< 2 (+ nx ny)) (pos? nx) (pos? ny))
                 (and (not equal-variances?)
                      (> nx 1) (> ny 1))) "Not enough observations.")
     (when paired? (assert (== nx ny) "Lengths of xs and ys should be equal."))
     (if paired?
       (-> (t-test-one-sample (map m/fast- xs ys) params)
           (assoc :paired? true))
       (let [{:keys [test-type ^double stat ^double alpha ^double df ^double mu ^double stderr]
              :as res} (test-two-samples-not-paired xs ys params)
             pvals (-> (test-pvalue-ci (r/distribution :t {:degrees-of-freedom df}) test-type stat alpha)
                       (update :confidence-interval (partial test-update-ci mu stderr)))]
         (assoc (merge pvals res) :t stat))))))

(def ^{:deprecated "Use [[t-test-two-samples]]"} ttest-two-samples t-test-two-samples)

(defn z-test-two-samples
  "Two samples z-test

  * `alpha` - significance level (default: `0.05`)
  * `sides` - one of: `:two-sided` (default), `:one-sided-less` (short: `:one-sided`) or `:one-sided-greater`
  * `mu` - mean (default: `0.0`)
  * `paired?` - unpaired or paired test, boolean (default: `false`)
  * `equal-variances?` - unequal or equal variances, boolean (default: `false`)"
  {:metadoc/categories #{:test}}
  ([xs ys] (z-test-two-samples xs ys {}))
  ([xs ys {:keys [paired? equal-variances?]
           :or {paired? false equal-variances? false}
           :as params}]
   (let [nx (count xs)
         ny (count ys)]
     (assert (or (and equal-variances? (< 2 (+ nx ny)) (pos? nx) (pos? ny))
                 (and (not equal-variances?)
                      (> nx 1) (> ny 1))) "Not enough observations.")
     (when paired? (assert (== nx ny) "Lengths of xs and ys should be equal."))
     (if paired?
       (-> (z-test-one-sample (map m/fast- xs ys) params)
           (assoc :paired? true))
       (let [{:keys [test-type ^double stat ^double alpha ^double ^double mu ^double stderr]
              :as res} (test-two-samples-not-paired xs ys params)
             pvals (-> (test-pvalue-ci r/default-normal test-type stat alpha)
                       (update :confidence-interval (partial test-update-ci mu stderr)))]
         (-> (merge res pvals)
             (assoc :z stat)
             (dissoc :df)))))))

(defn f-test
  "Variance F-test of two samples.

  * `alpha` - significance level (default: `0.05`)
  * `sides` - one of: `:two-sided` (default), `:one-sided-less` (short: `:one-sided`) or `:one-sided-greater` "
  {:metadoc/categories #{:test}}
  ([xs ys] (f-test xs ys {}))
  ([xs ys {:keys [sides ^double alpha]
           :or {sides :two-sided alpha 0.05}}]
   (let [nx (count xs)
         ny (count ys)
         dfx (dec nx)
         dfy (dec ny)
         F (/ (variance xs) (variance ys))
         distr (r/distribution :f {:denominator-degrees-of-freedom dfy
                                   :numerator-degrees-of-freedom dfx})]
     {:F F
      :stat F
      :estimate F
      :df [dfx dfy]
      :n [nx ny] :nx nx :ny ny
      :sides sides
      :test-type sides
      :p-value (p-value distr F sides)
      :confidence-interval (sides-case sides
                                       [(/ F ^double (r/icdf distr (- 1.0 (* alpha 0.5))))
                                        (/ F ^double (r/icdf distr (* alpha 0.5)))]
                                       [(/ F ^double (r/icdf distr (- 1.0 alpha))) ##Inf]
                                       [0.0 (/ F ^double (r/icdf distr alpha))])})))

(defn- pdt-gof
  [xs p ^double lambda]
  (let [cnt (count xs)
        n (double (reduce m/fast+ xs))
        df (dec cnt)
        p (or p (repeat cnt 1.0))
        psum (sum p)
        p (map (fn [^double p] (/ p psum)) p)
        xhat (map (fn [^double p] (* n p)) p)
        stat (condp = (double lambda)
               0.0 (* 2.0 ^double (reduce m/fast+ (map (fn [^long a ^double b]
                                                         (* a (- (m/log a) (m/log b)))) xs xhat)))
               -1.0 (* 2.0 ^double (reduce m/fast+ (map (fn [^double a ^long b]
                                                          (* a (- (m/log a) (m/log b)))) xhat xs)))
               (* (/ 2.0 (* lambda (inc lambda)))
                  ^double (reduce m/fast+ (map (fn [^long a ^double b]
                                                 (* a (dec (m/pow (/ a b) lambda)))) xs xhat))))]
    {:stat stat :df df :n n :expected xhat :p p
     :estimate (map (fn [^double v] (/ v n)) xs)}))

(defn- pdt-bootstrap-ci
  [{:keys [estimate ^long n ^double alpha ci-sides]} samples]
  (let [alpha (if (#{:both :two-sided} ci-sides) alpha (* alpha 2.0))
        d (r/distribution :multinomial {:trials n :ps (if (map? estimate) (vals estimate) estimate)})
        rands (apply map (comp m/seq->double-array vector) (r/->seq d samples))
        vs (sides-case ci-sides
                       (let [qs [(/ alpha 2.0) (- 1.0 (/ alpha 2.0))]]
                         (map #(v/div (quantiles % qs) (double n)) rands))
                       (let [q (/ alpha 2.0)]
                         (map #(vector (/ (quantile % q) n) 1.0) rands))
                       (let [q (- 1.0 (/ alpha 2.0))]
                         (map #(vector 0.0 (/ (quantile % q) n)) rands)))]
    (if (map? estimate) (zipmap (keys estimate) vs) vs)))

(defn- pdt-multi
  [ct ^double lambda]
  (let [xs (infer-ct ct)
        {:keys [^long n cols rows]} (contingency-table->marginals xs)
        n (double n)
        xhat (->> (for [[k1 ^long v1] rows
                        [k2 ^long v2] cols]
                    [[k1 k2] (/ (* v1 v2) n)])
                  (into {}))
        n1 (count (map first rows))
        n2 (count (map second cols))
        df (* (dec n1) (dec n2))
        stat (condp = (double lambda)
               0.0 (* 2.0 ^double (reduce (fn [^double sum [k ^long cnt]]
                                            (+ sum (* cnt (- (m/log cnt) (m/log (xhat k)))))) 0.0 xs))
               -1.0 (* 2.0 ^double (reduce (fn [^double sum [k ^double xhv]]
                                             (let [^double cnt (get xs k 0.0)]
                                               (+ sum (* xhv (- (m/log xhv) (m/log cnt) ))))) 0.0 xhat))
               (* (/ 2.0 (* lambda (inc lambda)))
                  ^double (reduce (fn [^double sum [k ^long cnt]]
                                    (+ sum (* cnt (dec (m/pow (/ cnt ^double (xhat k)) lambda))))) 0.0 xs)))]
    {:stat stat :df df :n n :k n1 :r n2
     :expected xhat
     :estimate (into {} (map (fn [[k ^long v]] [k (/ v n)]) xs))}))

(defn power-divergence-test
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {}))
  ([contingency-table-or-xs {:keys [^double lambda ci-sides sides p ^double alpha ^long bootstrap-samples
                                    ^long ddof]
                             :or {lambda m/TWO_THIRD sides :one-sided-greater ci-sides :two-sided
                                  alpha 0.05 bootstrap-samples 1000 ddof 0}}]
   (let [{:keys [df stat] :as res} (-> (if (and (sequential? contingency-table-or-xs)
                                                (every? number? contingency-table-or-xs))
                                         (pdt-gof contingency-table-or-xs p lambda)
                                         (pdt-multi contingency-table-or-xs lambda))
                                       (update :df (fn [^long df] (- df ddof))))
         distr (r/distribution :chi-squared {:degrees-of-freedom df})
         res (assoc res :lambda lambda :sides sides :test-type sides :ci-sides ci-sides :chi2 stat :alpha alpha :level (- 1.0 alpha)
                    :p-value (p-value distr stat sides))]
     (assoc res :confidence-interval (pdt-bootstrap-ci res bootstrap-samples)))))

(defn chisq-test
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda 1.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda 1.0))))

(defn multinomial-likelihood-ratio-test
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda 0.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda 0.0))))

(defn minimum-discrimination-information-test
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda -1.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda -1.0))))

(defn neyman-modified-chisq-test
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda -2.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda -2.0))))

(defn freeman-tukey-test
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda -0.5}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda -0.5))))

(defn cressie-read-test
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda m/TWO_THIRD}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda m/TWO_THIRD))))

;;

(defn- anova
  [xss]
  (let [Ni (map count xss)
        Zi (map mean xss)
        Z (/ (sum (mapcat identity xss)) (sum Ni))
        SSt (sum (map (fn [^double n ^double z]
                        (* n (m/sq (- z Z)))) Ni Zi))
        SSe (sum (map (fn [xs ^double zi]
                        (v/magsq (map (fn [^double v]
                                        (- v zi)) xs))) xss Zi))
        k (count Ni)
        DFt (dec k)
        DFe (- (sum Ni) k)
        MSe (/ SSe DFe)]
    {:n Ni :SSt SSt :SSe SSe :DFt DFt :DFe (int DFe) :MSt (/ SSt DFt) :MSe MSe}))

(defn- update-f-p-value
  [{:keys [DFt DFe ^double MSt ^double MSe] :as aov} sides]
  (let [F (/ MSt MSe)
        distr (r/distribution :f {:numerator-degrees-of-freedom DFt
                                  :denominator-degrees-of-freedom DFe})]
    (assoc aov
           :F F :df [DFt DFe] :stat F
           :p-value (p-value distr F sides))))

(defn one-way-anova-test
  ([xss] (one-way-anova-test xss {}))
  ([xss {:keys [sides]
         :or {sides :one-sided-greater}}]
   (update-f-p-value (anova xss) sides)))

(defn levene-test
  ([xss] (levene-test xss {}))
  ([xss {:keys [sides statistic scorediff]
         :or {sides :one-sided-greater statistic mean scorediff abs}}]
   (let [res (update-f-p-value (anova (map (fn [xs]
                                             (let [^double s (statistic xs)]
                                               (map (fn [^double v]
                                                      (scorediff (- v s))) xs))) xss)) sides)]
     (-> (assoc res :W (:F res))
         (dissoc :F)))))

(defn brown-forsythe-test
  ([xss] (levene-test xss {:statistic median}))
  ([xss params] (levene-test xss (assoc params :statistic median))))

(defn fligner-killeen-test
  ([xss] (fligner-killeen-test xss {}))
  ([xss {:keys [sides]
         :or {sides :one-sided-greater}}]
   (let [Z (mapcat (fn [xs]
                     (let [s (median xs)]
                       (map (fn [^double v] (abs (- v s))) xs))) xss)
         ranks (m/rank Z)
         rden (/ (* 2.0 (inc (count ranks))))
         qij (map (fn [^double r]
                    (r/icdf r/default-normal (+ 0.5 (* rden (inc r))))) ranks)
         {:keys [^double SSt ^double SSe ^int DFt ^int DFe]
          :as res} (->> (map count xss)
                        (reductions clojure.core/+ 0)
                        (partition 2 1)
                        (map (fn [[d t]] (drop d (take t qij))))
                        (anova))
         y (/ SSt SSe)
         chi2 (/ (* y (+ DFt DFe)) (inc y))
         distr (r/distribution :chi-squared {:degrees-of-freedom DFt})]
     (assoc res
            :chi2 chi2 :df DFt :stat chi2
            :p-value (p-value distr chi2 sides)))))

;; ad/ks tests

(defn- a2-stat
  ^double [^doubles xs d]
  (let [n (alength xs)]
    (reduce m/fast- (- n) (map (fn [^long idx]
                                 (* (/ (+ idx idx 1.0) n)
                                    (+ (m/log (r/cdf d (aget xs idx)))
                                       (m/log (r/ccdf d (aget xs (- n idx 1))))))) (range n)))))

(defn ad-test-one-sample
  "Anderson-Darling test"
  ([xs] (ad-test-one-sample xs r/default-normal))
  ([xs distribution-or-ys] (ad-test-one-sample xs distribution-or-ys {}))
  ([xs distribution-or-ys {:keys [sides kernel bandwidth] :or {sides :one-sided-greater kernel :gaussian}}]
   (let [d (cond
             (r/distribution? distribution-or-ys) distribution-or-ys
             (= kernel :enumerated) (r/distribution :enumerated-real {:data distribution-or-ys})
             :else (r/distribution :continuous-distribution {:data distribution-or-ys :kde kernel
                                                             :bandwidth bandwidth}))
         axs (m/seq->double-array (sort xs))
         stat (a2-stat axs d)
         n (alength axs)
         distr (r/distribution :anderson-darling {:n n})]
     {:stat stat :A2 stat :mean (mean axs) :stddev (stddev axs) :n n :sides sides
      :p-value (p-value distr stat sides)})))

(defn ks-test-one-sample
  "One sample Kolmogorov-Smirnov test"
  ([xs] (ks-test-one-sample xs r/default-normal))
  ([xs distribution-or-ys] (ks-test-one-sample xs distribution-or-ys {}))
  ([xs distribution-or-ys {:keys [sides kernel bandwidth distinct?]
                           :or {sides :two-sided kernel :gaussian distinct? true}}]
   (let [d (cond
             (r/distribution? distribution-or-ys) distribution-or-ys
             (= kernel :enumerated) (r/distribution :enumerated-real {:data distribution-or-ys})
             :else (r/distribution :continuous-distribution {:data distribution-or-ys :kde kernel
                                                             :bandwidth bandwidth}))
         xs (if distinct? (distinct xs) xs)
         n (count xs)
         dn (/ (double n))
         idxs (map (fn [^long i] (* i dn)) (range (inc n)))
         cdfs (map (partial r/cdf d) (sort xs))
         ^double dp (reduce m/fast-max (map m/fast- (rest idxs) cdfs))
         dn (- ^double (reduce m/fast-min (map m/fast- (butlast idxs) cdfs)))
         d (max dp dn)]
     {:n n :dp dp :dn dn :d d :sides sides
      :stat (sides-case sides d dp dn) 
      :p-value (sides-case sides
                           (p-value (r/distribution :kolmogorov-smirnov {:n n}) d :right)
                           (p-value (r/distribution :kolmogorov-smirnov+ {:n n}) dp :right)
                           (p-value (r/distribution :kolmogorov-smirnov+ {:n n}) dn :right))})))

(defn ks-test-two-samples
  "Two samples Kolmogorov-Smirnov test"
  ([xs ys] (ks-test-two-samples xs ys {}))
  ([xs ys {:keys [sides distinct?] :or {sides :two-sided distinct? true}}]
   (let [xs (if distinct? (distinct xs) xs)
         ys (if distinct? (distinct ys) ys)
         nx (count xs)
         ny (count ys)
         sort-idxs (m/order (concat xs ys))
         pdf-diffs (map (vec (concat (repeat nx (/ 1.0 nx)) (repeat ny (/ -1.0 ny)))) sort-idxs)
         [^double dn dp] (extent (reductions m/fast+ pdf-diffs))
         dn (- dn)
         dp (double dp) ;; 
         d (max dn dp)
         n (/ (* nx ny) (double (+ nx ny)))
         stat (* (m/sqrt n) (sides-case sides d dp dn))]
     {:n n :nx nx :ny ny :dp dp :dn dn :d d
      :stat stat
      :KS stat
      :sides sides
      :p-value (sides-case sides
                           (p-value (r/distribution :kolmogorov) stat :right)
                           (m/exp (* -2.0 (* stat stat)))
                           (m/exp (* -2.0 (* stat stat))))})))

(defn kruskal-test
  "Kruskal-Wallis rank sum test."
  ([xss] (kruskal-test xss {}))
  ([xss {:keys [sides] :or {sides :right}}] ;; as in R
   (let [k (count xss)
         df (dec k)
         xs (flatten xss)
         groups (mapcat (fn [[xs id]] (repeat (count xs) id)) (map vector xss (range)))
         n (count xs)
         r (m/rank1 xs)
         ties (vals (frequencies xs))
         ^double stat (->> (map vector r groups)
                           (group-by second)
                           (vals)
                           (map #(let [ranks (map first %)]
                                   (/ (m/sq (sum ranks)) (count ranks))))
                           (reduce m/fast+))
         stat (/ (- (/ (* 12.0 stat)
                       (* n (inc n)))
                    (* 3.0 (inc n)))
                 (- 1.0 (/ (sum (map (fn [^long t] (- (m/cb t) t)) ties))
                           (- (m/cb n) n))))]
     {:stat stat :n n :df df :k k :sides sides
      :p-value (p-value (r/distribution :chi-squared {:degrees-of-freedom df}) stat sides)})))

#_(m/unuse-primitive-operators)
