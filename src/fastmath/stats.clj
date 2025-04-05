(ns fastmath.stats
  "#### Statistics functions.

  * Descriptive statistics.
  * Correlation / covariance
  * Outliers
  * Confidence intervals
  * Extents
  * Effect size
  * Tests
  * Histogram
  * ACF/PACF
  * Bootstrap (see `fastmath.stats.bootstrap`)
  * Binary measures

  Functions are backed by Apache Commons Math or SMILE libraries. All work with Clojure sequences.

  ##### Descriptive statistics

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
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.distance :as d]
            [fastmath.vector :as v]
            [fastmath.interpolation.step :as step-interp]
            [fastmath.interpolation.linear :as linear-interp]
            [fastmath.optimization.lbfgsb :as lbfgsb]
            [fastmath.optimization :as opt]
            [fastmath.kernel.density :as kd]
            [fastmath.special :as special]
            [fastmath.solver :as solver])
  (:import [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.stat.descriptive.rank Percentile Percentile$EstimationType]
           [org.apache.commons.math3.stat.descriptive.moment Kurtosis Skewness]
           [org.apache.commons.math3.stat.correlation KendallsCorrelation SpearmansCorrelation PearsonsCorrelation]
           [org.apache.commons.math3.stat.regression SimpleRegression]
           [org.apache.commons.math3.analysis.integration RombergIntegrator]
           [org.apache.commons.math3.analysis UnivariateFunction]
           [fastmath.java Array]
           [fastmath.vector Vec2 Vec3]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators)

(def ^{:doc "List of estimation strategies for [[percentile]]/[[quantile]] functions."}
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

(defn- kahan-step
  ^Vec2 [^Vec2 b ^double v]
  (let [av (m/- v (.y b))
        nsum (m/+ av (.x b))]
    (Vec2. nsum (m/- nsum (.x b) av))))

(defn- neumayer-step
  ^Vec2 [^Vec2 b ^double v]
  (let [t (m/+ (.x b) v)]
    (Vec2. t (m/+ (.y b) (if (m/>= (m/abs (.x b)) (m/abs v))
                           (m/+ (m/- (.x b) t) v)
                           (m/+ (m/- v t) (.x b)))))))

(defn- klein-step
  ^Vec3 [^Vec3 b ^double v]
  (let [t (m/+ (.x b) v)
        c (if (m/>= (m/abs (.x b)) (m/abs v))
            (m/+ (m/- (.x b) t) v)
            (m/+ (m/- v t) (.x b)))
        nt (m/+ (.y b) c)
        cc (if (m/>= (m/abs (.y b)) (m/abs c))
             (m/+ (m/- (.y b) nt) c)
             (m/+ (m/- c nt) (.y b)))]
    (Vec3. t nt (m/+ (.z b) cc))))

(defn sum
  "Sum of all `vs` values.

   Possible compensated summation methods are: `:kahan`, `:neumayer` and `:klein`"
  (^double [vs]
   (if (= (type vs) m/double-array-type)
     (Array/sum ^doubles vs)
     (reduce + vs)))
  (^double [vs compensation-method]
   (case compensation-method
     :kahan (let [^Vec2 r (reduce kahan-step (Vec2. 0.0 0.0) vs)] (.x r))
     :neumayer (v/sum (reduce neumayer-step (Vec2. 0.0 0.0) vs))
     :klein (v/sum (reduce klein-step (Vec3. 0.0 0.0 0.0) vs))
     (sum vs))))

(defn percentile
  "Calculate percentile of a `vs`.

  Percentile `p` is from range 0-100.

  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html).

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values: `:legacy`, `:r1` to `:r9`. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html)

  See also [[quantile]]."
  (^double [vs ^double p]
   (StatUtils/percentile (m/seq->double-array vs) p))
  (^double [vs ^double p estimation-strategy]
   (let [^Percentile perc (.withEstimationType (Percentile.) (get estimation-strategies-list estimation-strategy Percentile$EstimationType/LEGACY))]
     (.evaluate perc (m/seq->double-array vs) p))))

(defn percentiles
  "Calculate percentiles of a `vs`.

  Percentiles are sequence of values from range 0-100.

  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html).

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values: `:legacy`, `:r1` to `:r9`. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html)

  See also [[quantile]]."
  ([vs] (percentiles vs [25 50 75 100]))
  ([vs ps] (percentiles vs ps nil))
  ([vs ps estimation-strategy]
   (let [^Percentile perc (.withEstimationType (Percentile.) (or (estimation-strategies-list estimation-strategy) Percentile$EstimationType/LEGACY))]
     (.setData perc (m/seq->double-array vs))
     (mapv (fn [^double p] (.evaluate perc p)) ps))))

(defn quantile
  "Calculate quantile of a `vs`.

  Quantile `q` is from range 0.0-1.0.

  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html) for interpolation strategy.

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values: `:legacy`, `:r1` to `:r9`. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html)

  See also [[percentile]]."
  (^double [vs ^double q]
   (percentile vs (m/constrain (* q 100.0) 0.0 100.0)))
  (^double [vs ^double q estimation-strategy]
   (percentile vs (m/constrain (* q 100.0) 0.0 100.0) estimation-strategy)))

(defn quantiles
  "Calculate quantiles of a `vs`.

  Quantilizes is sequence with values from range 0.0-1.0.

  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html) for interpolation strategy.

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values: `:legacy`, `:r1` to `:r9`. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html)

  See also [[percentiles]]."
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
        wsum (sum probabilities)
        weights (conj (reductions + (map (fn [^double p] (/ p wsum)) probabilities)) 0.0)]
    (case method
      :linear (linear-interp/linear weights data)
      :average (let [interp1 (step-interp/step-before weights data)
                     interp2 (step-interp/step-after weights data)]
                 (fn [^double x] (* 0.5 (+ ^double (interp1 x)
                                          ^double (interp2 x)))))
      :step (step-interp/step-before weights data))))

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
  (^double [vs estimation-strategy]
   (percentile vs 50.0 estimation-strategy))
  (^double [vs]
   (percentile vs 50.0)))

(defn median-3
  "Median of three values. See [[median]]."
  ^double [^double a ^double b ^double c]
  (m/max (m/min a b) (m/min (m/max a b) c)))

(defn mean
  "Calculate mean of `vs` with optional `weights`."
  (^double [vs] (StatUtils/mean (m/seq->double-array vs)))
  (^double [vs weights] (/ (sum (map * vs weights)) (sum weights))))

(defn geomean
  "Geometric mean for positive values only with optional `weights`"
  (^double [vs] (m/exp (mean (map (fn [^double v] (m/log v)) vs))))
  (^double [vs weights] (m/exp (mean (map (fn [^double v] (m/log v)) vs) weights))))

(defn harmean
  "Harmonic mean with optional `weights`"
  (^double [vs] (/ (mean (map (fn [^double v] (/ v)) vs))))
  (^double [vs weights] (/ (mean (map (fn [^double v] (/ v)) vs) weights))))

(defn powmean
  "Generalized power mean"
  (^double [vs ^double power]
   (cond
     (zero? power) (geomean vs)
     (m/one? power) (mean vs)
     (== power m/THIRD) (m/cb (mean (map #(m/cbrt %) vs)))
     (== power 0.5) (m/sq (mean (map #(m/sqrt %) vs)))
     (== power 2.0) (m/sqrt (mean (map m/sq vs)))
     (== power 3.0) (m/cbrt (mean (map m/cb vs))) 
     :else (m/pow (mean (map (fn [^double v] (m/pow v power)) vs)) (/ power))))
  (^double [vs weights ^double power]
   (cond
     (zero? power) (geomean vs weights)
     (m/one? power) (mean vs weights)
     (== power m/THIRD) (m/cb (mean (map m/cbrt vs) weights))
     (== power 0.5) (m/sq (mean (map m/sqrt vs) weights))
     (== power 2.0) (m/sqrt (mean (map m/sq vs) weights))
     (== power 3.0) (m/cbrt (mean (map m/cb vs) weights)) 
     :else (m/pow (mean (map (fn [^double v] (m/pow v power)) vs) weights) (/ power)))))

(defn wmean
  "Weighted mean"
  {:deprecated "Use `mean`"}
  (^double [vs] (mean vs))
  (^double [vs weights]
   (/ (sum (map * vs weights)) (sum weights))))

(defn population-variance
  "Calculate population variance of `vs`.

  See [[variance]]."
  (^double [vs]
   (StatUtils/populationVariance (m/seq->double-array vs)))
  (^double [vs ^double mu]
   (StatUtils/populationVariance (m/seq->double-array vs) mu)))

(defn population-wvariance
  "Calculate weighted population variance of `vs`."
  ^double [vs freqs]
  (let [sw (sum freqs)
        mu (/ (sum (map * vs freqs)) sw)
        v (sum (map (fn [^double x ^double w]
                      (* w (m/sq (- x mu)))) vs freqs))]
    (/ v sw)))

(defn variance
  "Calculate variance of `vs`.

  See [[population-variance]]."
  (^double [vs]
   (StatUtils/variance (m/seq->double-array vs)))
  (^double [vs ^double mu]
   (StatUtils/variance (m/seq->double-array vs) mu)))

(defn wvariance
  "Calculate weighted (unbiased) variance of `vs`."
  ^double [vs freqs]
  (let [sw (sum freqs)
        mu (/ (sum (map * vs freqs)) sw)
        v (sum (map (fn [^double x ^double w]
                      (* w (m/sq (- x mu)))) vs freqs))]
    (/ v (dec sw))))

(defn population-stddev
  "Calculate population standard deviation of `vs`.

  See [[stddev]]."
  (^double [vs]
   (m/sqrt (population-variance vs)))
  (^double [vs ^double mu]
   (m/sqrt (population-variance vs mu))))

(defn population-wstddev
  "Calculate population weighted standard deviation of `vs`"
  ^doubles [vs weights]
  (m/sqrt (population-wvariance vs weights)))

(defn stddev
  "Calculate standard deviation of `vs`.

  See [[population-stddev]]."
  (^double [vs]
   (m/sqrt (variance vs)))
  (^double [vs ^double mu]
   (m/sqrt (variance vs mu))))

(defn wstddev
  "Calculate weighted (unbiased) standard deviation of `vs`"
  ^doubles [vs freqs]
  (m/sqrt (wvariance vs freqs)))

(defn variation
  "Coefficient of variation CV = stddev / mean"
  ^double [vs]
  (let [vs (m/seq->double-array vs)]
    (/ (stddev vs)
       (mean vs))))

(defn median-absolute-deviation
  "Calculate MAD"
  (^double [vs] (median-absolute-deviation vs nil))
  (^double [vs center]
   (let [m (double (or center (median vs)))]
     (median (map (fn [^double x] (m/abs (- x m))) vs))))
  (^double [vs center estimation-strategy]
   (let [m (double (or center (median vs estimation-strategy)))]
     (median (map (fn [^double x] (m/abs (- x m))) vs) estimation-strategy))))

(def ^{:doc "Alias for [[median-absolute-deviation]]"}
  mad median-absolute-deviation)

(defn mean-absolute-deviation
  "Calculate mean absolute deviation"
  (^double [vs] (mean-absolute-deviation vs nil))
  (^double [vs center]
   (let [m (double (or center (mean vs)))]
     (mean (map (fn [^double x] (m/abs (- x m))) vs)))))

(defn sem
  "Standard error of mean"
  ^double [vs]
  (let [s (stddev vs)]
    (/ s (m/sqrt (count vs)))))

(defmacro ^:private build-extent
  [nm mid ext]
  `(defn ~nm
     ~(str " -/+ " ext " and " mid)
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

  Returns a sequence of outliers.

  Optional `estimation-strategy` argument can be set to change quantile calculations estimation type. See [[estimation-strategies]]."
  ([vs]
   (outliers vs :legacy))
  ([vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         q1 (percentile avs 25.0 estimation-strategy)
         q3 (percentile avs 75.0 estimation-strategy)]
     (outliers avs q1 q3)))
  ([vs ^double q1 ^double q3]
   (let [iqr (* 1.5 (- q3 q1))
         lif-thr (- q1 iqr)
         uif-thr (+ q3 iqr)]
     (filter (fn [^double v]
               (or (< v lif-thr)
                   (> v uif-thr))) vs))))

(defn remove-outliers
  "Remove outliers defined as values outside inner fences.

  Let Q1 is 25-percentile and Q3 is 75-percentile. IQR is `(- Q3 Q1)`.

  * LIF (Lower Inner Fence) equals `(- Q1 (* 1.5 IQR))`.
  * UIF (Upper Inner Fence) equals `(+ Q3 (* 1.5 IQR))`.

  Returns a sequence without outliers.

  Optional `estimation-strategy` argument can be set to change quantile calculations estimation type. See [[estimation-strategies]]."
  ([vs]
   (outliers vs :legacy))
  ([vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         q1 (percentile avs 25.0 estimation-strategy)
         q3 (percentile avs 75.0 estimation-strategy)]
     (outliers avs q1 q3)))
  ([vs ^double q1 ^double q3]
   (let [iqr (* 1.5 (- q3 q1))
         lif-thr (- q1 iqr)
         uif-thr (+ q3 iqr)]
     (remove (fn [^double v]
               (or (< v lif-thr)
                   (> v uif-thr))) vs))))

(defn wmodes
  "Returns weighted modes.

  `vs` can contain anything, in case of ties returns all modes."
  ([vs] (wmodes vs (repeat 1.0)))
  ([vs weights]
   (let [all (->> (map vector vs weights)
                  (reduce (fn [b [v w]]
                            (update b v (fnil m/+ 0.0) w)) {})
                  (sort-by second m/>))
         best-score (-> all first second double)]
     (->> all
          (take-while (fn [p] (m/== best-score (double (second p)))))
          (map first)))))

(defn wmode
  "Returns weighted mode.

  `vs` can contain anything, in case of ties returns only one (possibly random) mode."
  ([vs] (wmode vs (repeat 1.0)))
  ([vs weights]
   (first (wmodes vs weights))))

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
       :kde (let [kde (kd/kernel-density (get opts :kernel :gaussian) avs opts)]
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
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (Array/min ^doubles vs)
    (reduce m/min vs)))

(defn maximum
  "Maximum value from sequence."
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (Array/max ^doubles vs)
    (reduce m/max vs)))

(defn span
  "Width of the sample, maximum value minus minimum value"
  ^double [vs]
  (let [avs (m/seq->double-array vs)]
    (- (maximum avs) (minimum avs))))

(defn extent
  "Return extent (min, max, mean) values from sequence. Mean is optional (default: true)"
  ([vs] (extent vs true))
  ([vs mean?]
   (let [^double fv (first vs)
         mm (reduce (fn [^Vec2 curr ^double v]
                      (Vec2. (min (.x curr) v) (max (.y curr) v))) (Vec2. fv fv) (rest vs))]
     (if mean? (conj mm (mean vs)) mm))))

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
             (m/one? order) m/identity-double
             (== order 2.0) m/sq
             (== order 3.0) m/cb
             (== order 4.0) (fn ^double [^double diff] (m/sq (m/sq diff)))
             :else (fn ^double [^double diff] (m/pow diff order)))
         a (if absolute? m/abs m/identity-double)]
     (loop [idx (int 0)]
       (when (< idx cin)
         (aset out idx ^double (f (a (- (aget in idx) center))))
         (recur (inc idx))))
     (/ (if mean? (mean out) (sum out)) nf))))

(def ^{:deprecated "Use [[moment]] function"} second-moment moment)

;;

(defn l-moment
  "Calculates L-moment, TL-moment (trimmed) or (T)L-moment ratios.

  Options:
  - `:s` (default: 0) - number of left trimmed values
  - `:r` (default: 0) - number of right tirmmed values
  - `:sorted?` (default: false) - if input is already sorted
  - `:ratio?` (default: false) - normalized l-moment, l-moment ratio"
  (^double [vs ^long order] (l-moment vs order nil))
  (^double [vs ^long order {:keys [^long s ^long t sorted? ratio?]
                            :or {s 0 t 0}
                            :as opts}]
   (if (m/zero? order)
     1.0
     (let [^doubles svs (m/seq->double-array (if sorted? vs (sort vs)))]
       (if ratio?
         (let [nopts (assoc opts :sorted? true :ratio? false)
               l2 (l-moment svs 2 nopts)]
           (m// (l-moment svs order nopts) l2))
         (let [r- (m/dec order)
               s+ (m/inc s)
               n (alength svs)
               n-t+ (m/- n t -1)]
           (m// (double (reduce (fn [^double b1 ^long k]
                                  (let [c1 (m/- (m/+ r- s) k)
                                        c2 (m/+ t k)]
                                    (-> (if (m/even? k) (m/combinations r- k) (m/- (m/combinations r- k)))
                                        (m/* (double (reduce (fn [^double b2 ^long j]
                                                               (let [j- (m/dec j)]
                                                                 (-> (m/* (m/combinations j- c1)
                                                                          (m/combinations (m/- n j) c2)
                                                                          (Array/aget svs j-))
                                                                     (m/+ b2)))) 0.0 (range s+ n-t+))))
                                        (m/+ b1))))
                                0.0 (range order)))
                (m/* order (m/combinations n (m/+ order s t))))))))))

(defn l-variation
  "Coefficient of L-variation, L-CV"
  ^double [vs]
  (let [svs (m/seq->double-array vs)]
    (java.util.Arrays/sort svs)
    (m// (l-moment svs 2 {:sorted? true})
         (mean svs))))

;;; expectile

(defn- expectile-map
  [vs ^double tau ^double tau- ^double t]
  (map (fn [^double v]
         (m/* (if (m/<= v t) tau- tau) (m/- t v))) vs))

(defn- expectile-target
  ([vs weights ^double tau ^double tau-]
   (fn ^double [^double t]
     (let [m (expectile-map vs tau tau- t)]
       (if weights (mean m weights) (mean m)))))  )

(defn expectile
  "Calculate expectile for given tau.

  If tau=0.5, returns mean."
  (^double [vs ^double tau] (expectile vs nil tau))
  (^double [vs weights ^double tau]
   (let [avg (if weights (mean vs weights) (mean vs))]
     (if (m/== tau 0.5)
       avg
       (let [[^double x0 ^double x1] (if (m/> tau 0.5)
                                       [avg (maximum vs)]
                                       [(minimum vs) avg])]
         (if (m/== x0 x1)
           x0
           (solver/find-root (expectile-target vs weights tau (m/- 1.0 tau)) x0 x1)))))))

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
     (->> vs
          (filter (fn [^double v]
                    (or (m/nan? v)
                        (<= low v high))))
          (map (fn [^double v] (if (m/nan? v) nan v)))))))

(defn trim-lower
  "Trim data below given quanitle, default: 0.2."
  ([vs] (trim-lower vs 0.2))
  ([vs quantile] (trim-lower vs quantile :legacy))
  ([vs ^double quantile estimation-strategy]
   (let [[q qmid] (quantiles (remove m/nan? vs) [quantile 0.5] estimation-strategy)]
     (trim vs q ##Inf qmid))))

(defn trim-upper
  "Trim data above given quanitle, default: 0.2."
  ([vs] (trim-upper vs 0.2))
  ([vs quantile] (trim-upper vs quantile :legacy))
  ([vs ^double quantile estimation-strategy]
   (let [[q qmid] (quantiles (remove m/nan? vs) [quantile 0.5] estimation-strategy)]
     (trim vs ##-Inf q qmid))))

;; On More Robust Estimation of Skewness and Kurtosis: Simulation and Application to the S&P500 Index

(defn- yule-skewness
  ^double [vs ^double u]
  (let [[^double q1 ^double q2 ^double q3] (quantiles vs [u 0.5 (- 1.0 u)])]
    (/ (+ q3 (* -2.0 q2) q1)
       (- q3 q1))))

(defn- bowley-skewness
  ^double [vs]
  (let [[^double q1 ^double q3] (quantiles vs [0.25 0.75])]
    (/ (+ q3 q1) (- q3 q1))))

(defn- hogg-skewness
  ^double [vs]
  (let [m25 (mean (trim vs 0.25))
        u005 (mean (trim-lower vs 0.95))
        l005 (mean (trim-upper vs 0.05))]
    (/ (- u005 m25) (- m25 l005))))

(defn skewness
  "Calculate skewness from sequence.

  Possible types: `:G1` (default), `:g1` (`:pearson`), `:b1`, `:B1` (`:yule`), `:B3`, `:skew`, `:mode`, `:bowley`, `:hogg` or `:median`."
  (^double [vs] (skewness vs :G1))
  (^double [vs typ]
   (let [vs (m/seq->double-array vs)]
     (if (sequential? typ)
       (cond
         (= :mode (first typ)) (let [[_ method opts] typ]
                                 (/ (- (mean vs) (mode vs method opts)) (stddev vs)))
         (#{:B1 :yule} (first typ)) (yule-skewness vs (second typ)))
       (case typ
         :mode (/ (- (mean vs) (mode vs)) (stddev vs))
         :median (/ (* 3.0 (- (mean vs) (median vs))) (stddev vs))
         :bowley (bowley-skewness vs)
         :hogg (hogg-skewness vs)
         (:B1 :yule) (yule-skewness vs 0.25)
         :B3 (let [v (median vs)]
               (/ (- (mean vs) v)
                  (moment vs 1.0 {:absolute? true :center v})))
         (let [^Skewness k (Skewness.)
               n (alength vs)
               v (.evaluate k vs)]
           (cond
             (= :b1 typ) (* v (/ (* (- n 2.0) (dec n)) (* n n)))
             (#{:pearson :g1} typ) (* v (/ (- n 2.0) (m/sqrt (* n (dec n)))))
             (= :skew typ) (* v (/ (- n 2.0) (* n (m/sqrt (dec n))))) ;; artificial, to match BCa skew definition
             :else v)))))))

;; centered
(defn- moors-kurtosis
  ^double [vs]
  (let [[^double e1 ^double e2 ^double e3
         ^double e5 ^double e6 ^double e7] (quantiles vs [0.125 0.25 0.375 0.625 0.75 0.875])]
    (- (/ (+ (- e7 e5) (- e3 e1))
          (- e6 e2)) 1.23)))

;; centered
(defn- crow-kurtosis
  (^double [vs] (crow-kurtosis vs 0.025 0.25))
  (^double [vs ^double alpha ^double beta]
   (let [[^double a1 ^double a2 ^double b1 ^double b2] (quantiles vs [alpha (- 1.0 alpha)
                                                                      beta (- 1.0 beta)])]
     (- (/ (- a2 a1) (- b2 b1)) 2.91))))

;; centered
(defn- hogg-kurtosis
  (^double [vs] (hogg-kurtosis vs 0.05 0.5))
  (^double [vs ^double alpha ^double beta]
   (let [ua (mean (trim-lower vs (- 1.0 alpha)))
         ub (mean (trim-lower vs (- 1.0 beta)))
         la (mean (trim-upper vs alpha))
         lb (mean (trim-upper vs beta))]
     (- (/ (- ua la) (- ub lb)) 2.59))))

(defn kurtosis
  "Calculate kurtosis from sequence.

  Possible typs: `:G2` (default), `:g2` (or `:excess`), `:geary`, ,`:crow`, `:moors`, `:hogg` or `:kurt`."
  (^double [vs] (kurtosis vs nil))
  (^double [vs typ]
   (let [vs (m/seq->double-array vs)
         n (alength vs)]
     (if (sequential? typ)
       (condp = (first typ)
         :crow (apply crow-kurtosis vs (rest typ))
         :hogg (apply hogg-kurtosis vs (rest typ)))
       (condp = typ
         :geary (/ (mean-absolute-deviation vs)
                   (population-stddev vs))
         :moors (moors-kurtosis vs)
         :crow (crow-kurtosis vs)
         :hogg (hogg-kurtosis vs)
         (let [^Kurtosis k (Kurtosis.)
               v (.evaluate k vs)]
           (cond
             (#{:excess :g2} typ) (/ (- (/ (* v (- n 2) (- n 3)) (dec n)) 6.0)
                                     (inc n))
             (= :kurt typ) (+ 3.0 (/ (- (/ (* v (- n 2) (- n 3)) (dec n)) 6.0)
                                     (inc n)))
             :else v)))))))

(defn ci
  "T-student based confidence interval for given data. Alpha value defaults to 0.05.

  Last value is mean."
  ([vs] (ci vs 0.05))
  ([vs ^double alpha]
   (let [vsa (m/seq->double-array vs)
         cnt (count vs)
         dist (r/distribution :t {:degrees-of-freedom (dec cnt)})
         ^double crit-val (r/icdf dist (- 1.0 (* 0.5 alpha)))
         mean-ci (/ (* crit-val (stddev vsa)) (m/sqrt cnt))
         mn (mean vsa)]
     [(- mn mean-ci) (+ mn mean-ci) mn])))

;; https://ocw.mit.edu/courses/mathematics/18-05-introduction-to-probability-and-statistics-spring-2014/readings/MIT18_05S14_Reading24.pdf
(defn bootstrap-ci
  "Bootstrap method to calculate confidence interval.

  Alpha defaults to 0.98, samples to 1000.
  Last parameter is statistical function used to measure, default: [[mean]].

  Returns ci and statistical function value."
  {:deprecated "Please use fastmath.stats.boostrap/ci-basic instead"}
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
  {:doc "Generate set of samples of given size from provided data.

  Default `samples` is 200, number of `size` defaults to sample size."
   :deprecated "Please use fastmath.stats.bootstrap/bootstrap instead"}
  ([vs] (bootstrap vs 200))
  ([vs samples] (bootstrap vs samples (count vs)))
  ([vs samples size]
   (let [dist (r/distribution :enumerated-real {:data vs})]
     (repeatedly samples #(r/->seq dist size)))))

(defn stats-map
  "Calculate several statistics of `vs` and return as map.

  Optional `estimation-strategy` argument can be set to change quantile calculations estimation type. See [[estimation-strategies]]."
  ([vs] (stats-map vs :legacy))
  ([vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         sz (alength avs)
         mn (Array/min avs)
         mx (Array/max avs)
         sm (Array/sum avs)
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
  [vs]
  (let [m (mean vs)]
    (map (fn [^double v]
           (- v m)) vs)))

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
  (^double [[vs1 vs2]] (covariance vs1 vs2))
  (^double [vs1 vs2]
   (let [avs1 (m/seq->double-array vs1)
         avs2 (m/seq->double-array vs2)]
     (/ (v/dot (v/shift avs1 (- (mean avs1)))
               (v/shift avs2 (- (mean avs2))))
        (dec (alength avs1))))))

(defn correlation
  "Correlation of two sequences."
  (^double [[vs1 vs2]] (correlation vs1 vs2))
  (^double [vs1 vs2]
   (let [avs1 (m/seq->double-array vs1)
         avs2 (m/seq->double-array vs2)
         cov (covariance avs1 avs2)
         v1 (variance avs1)
         v2 (variance avs2)]
     (if (or (zero? v1) (zero? v2))
       ##NaN
       (/ cov (m/sqrt (* v1 v2)))))))

(defn spearman-correlation
  "Spearman's correlation of two sequences."
  (^double [[vs1 vs2]] (spearman-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation (SpearmansCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn pearson-correlation
  "Pearson's correlation of two sequences."
  (^double [[vs1 vs2]] (pearson-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation (PearsonsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn kendall-correlation
  "Kendall's correlation of two sequences."
  (^double [[vs1 vs2]] (kendall-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation (KendallsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn ^{:deprecated "Use [[dissimilarity]]."} kullback-leibler-divergence
  "Kullback-Leibler divergence of two sequences."
  (^double [[vs1 vs2]] (kullback-leibler-divergence vs1 vs2))
  (^double [vs1 vs2]
   (let [res (->> (map vector vs1 vs2)
                  (remove #(zero? (v/prod %)))
                  (map (fn [[^double p ^double q]] (* p (m/log (/ p q))))))]
     (if (seq res) (sum res) ##Inf))))

(defn ^{:deprecated "Use [[dissimilarity]]."} jensen-shannon-divergence
  "Jensen-Shannon divergence of two sequences."
  (^double [[vs1 vs2]] (jensen-shannon-divergence vs1 vs2))
  (^double [vs1 vs2]
   (let [m (v/mult (mapv + vs1 vs2) 0.5)]
     (* 0.5 (+ (kullback-leibler-divergence vs1 m)
               (kullback-leibler-divergence vs2 m))))))

(defn coefficient-matrix
  "Generate coefficient (correlation, covariance, any two arg function) matrix from seq of seqs. Row order.

  Default method: pearson-correlation"
  ([vss] (coefficient-matrix vss pearson-correlation))
  ([vss measure-fn] (coefficient-matrix vss measure-fn false))
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

  Possible measures: `:pearson` (default), `:kendall`, `:spearman`."
  ([vss] (correlation-matrix vss :pearson))
  ([vss measure]
   (let [measure (get {:pearson pearson-correlation
                       :kendall kendall-correlation
                       :spearman spearman-correlation} measure pearson-correlation)]
     (coefficient-matrix vss measure true))))

(defn covariance-matrix
  "Generate covariance matrix from seq of seqs. Row order."
  [vss] (coefficient-matrix vss covariance true))

(defn- maybe-number->seq [vs1 vs2] (if (number? vs2) [(seq vs1) (repeat (count vs1) vs2)] [vs1 vs2]))

(defn me
  "Mean error"
  (^double [[vs1 vs2-or-val]] (me vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (mean (map - v1 v2)))))

(defn mae
  "Mean absolute error"
  (^double [[vs1 vs2-or-val]] (mae vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (mean (map (comp m/abs -) v1 v2)))))

(defn mape
  "Mean absolute percentage error"
  (^double [[vs1 vs2-or-val]] (mape vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (mean (map (fn [^double a ^double b]
                  (m/abs (/ (- a b) a))) v1 v2)))))

(defn rss
  "Residual sum of squares"
  (^double [[vs1 vs2-or-val]] (rss vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (sum (map (comp m/sq -) v1 v2)))))

(defn r2
  "R2"
  (^double [[vs1 vs2-or-val]] (r2 vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (- 1.0 (/ (rss vs1 vs2-or-val)
             (moment vs1 2 {:mean? false}))))
  (^double [vs1 vs2-or-val ^double no-of-variables]
   (let [rr (r2 vs1 vs2-or-val)
         n (count vs1)]
     (- 1.0 (* (- 1.0 rr)
               (/ (dec n)
                  (- n no-of-variables 1.0)))))))

(defn mse
  "Mean squared error"
  (^double [[vs1 vs2-or-val]] (mse vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (mean (map (comp m/sq -) v1 v2)))))

(defn rmse
  "Root mean squared error"
  (^double [[vs1 vs2-or-val]] (rmse vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (m/sqrt (mse vs1 vs2-or-val))))

(defn count=
  "Count equal values in both seqs. Same as [[L0]]"
  (^long [[vs1 vs2-or-val]] (count= vs1 vs2-or-val))
  (^long [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (count (filter (fn [^double v] (zero? v)) (map - v1 v2))))))

(def ^{:doc "Count equal values in both seqs. Same as [[count==]]"} L0 count=)

(defn L1
  "Manhattan distance"
  (^double [[vs1 vs2-or-val]] (L1 vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (d/manhattan v1 v2))))

(defn L2sq
  "Squared euclidean distance"
  (^double [[vs1 vs2-or-val]] (L2sq vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (d/euclidean-sq v1 v2))))

(defn L2
  "Euclidean distance"
  (^double [[vs1 vs2-or-val]] (L2 vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (d/euclidean v1 v2))))

(defn LInf
  "Chebyshev distance"
  (^double [[vs1 vs2-or-val]] (LInf vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (d/chebyshev v1 v2))))

(defn psnr
  "Peak signal to noise, `max-value` is maximum possible value (default: max from `vs1` and `vs2`)"
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
  (^long [vs] (estimate-bins vs :freedman-diaconis))
  (^long [vs bins-or-estimate-method]
   (if-not (keyword? bins-or-estimate-method)
     (or bins-or-estimate-method (estimate-bins vs))
     (let [n (count vs)]
       (min n (int (case bins-or-estimate-method
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

(defn- constrain-data
  [vs ^double mn ^double mx]
  (filter (fn [^double v] (<= mn v mx)) vs))

(defn- process-vs-and-bins
  [vs bins-or-estimate-method ^double mn ^double mx]
  (if (sequential? bins-or-estimate-method)
    (let [[^double nmn ^double nmx] (extent bins-or-estimate-method)]
      [(constrain-data vs nmn nmx)
       nmn nmx (sort bins-or-estimate-method)])
    (let [nvs (constrain-data vs mn mx)
          bins (if (== mn mx) 1 (estimate-bins nvs bins-or-estimate-method))]
      [nvs mn mx (m/slice-range mn mx (inc bins))])))

(defn- histogram-internal
  [[vs ^double mn ^double mx intervals]]
  (let [diff (- mx mn)
        bins (if (zero? diff) 1 (dec (count intervals)))
        bins- (dec bins)
        step (/ diff bins)
        search-array (double-array intervals)
        buff (long-array bins)
        sum (double-array bins)
        size (count vs)
        dsize (double size)]

    (doseq [^double v vs]
      (let [b (java.util.Arrays/binarySearch ^doubles search-array v)
            pos (unchecked-int (min bins- (long (if (neg? b) (m/abs (+ b 2)) b))))]
        (fastmath.java.Array/inc ^longs buff pos)
        (fastmath.java.Array/add ^doubles sum pos v)))

    (let [bins-map (map (fn [[^double mn ^double mx] ^long cnt ^double s]
                          {:min mn :max mx :count cnt
                           :step (- mx mn)
                           :avg (/ s cnt)
                           :probability (/ cnt dsize)}) (partition 2 1 search-array) buff sum)]
      {:size bins
       :step step
       :samples size
       :min mn
       :max mx
       :bins (map (juxt :min :count) bins-map)
       :bins-maps bins-map
       :intervals intervals
       :frequencies (into {} (map (juxt :avg :count) bins-map))})))

(defn histogram
  "Calculate histogram.

  Estimation method can be a number, named method: `:sqrt` `:sturges` `:rice` `:doane` `:scott` `:freedman-diaconis` (default) or a sequence of points used as intervals.
  In the latter case or when `mn` and `mx` values are provided - data will be filtered to fit in desired interval(s).

  Returns map with keys:

  * `:size` - number of bins
  * `:step` - average distance between bins
  * `:bins` - seq of pairs of range lower value and number of elements
  * `:min` - min value
  * `:max` - max value
  * `:samples` - number of used samples
  * `:frequencies` - a map containing counts for bin's average
  * `:intervals` - intervals used to create bins
  * `:bins-maps` - seq of maps containing:
    * `:min` - lower bound
    * `:max` - upper bound
    * `:step` - actual distance between bins 
    * `:count` - number of elements
    * `:avg` - average value
    * `:probability` - probability for bin

  If difference between min and max values is `0`, number of bins is set to 1."
  ([vs] (histogram vs :freedman-diaconis))
  ([vs bins-or-estimate-method] (histogram vs bins-or-estimate-method (extent
                                                                       (if (sequential? (first vs))
                                                                         (flatten vs) vs))))
  ([vs bins-or-estimate-method [^double mn ^double mx]] (histogram vs bins-or-estimate-method mn mx))
  ([vs bins-or-estimate-method ^double mn ^double mx]
   (if (sequential? (first vs))
     (let [sbins? (sequential? bins-or-estimate-method)
           [^double nmn ^double nmx] (if sbins? (extent bins-or-estimate-method) [mn mx])
           nvs (map #(constrain-data % nmn nmx) vs)
           intervals (if sbins?
                       (sort bins-or-estimate-method)
                       (if (== nmn nmx)
                         [nmn nmx]
                         (->> nvs
                              (map #(estimate-bins % bins-or-estimate-method))
                              (reduce m/max 1.0)
                              (int)
                              (inc)
                              (m/slice-range nmn nmx))))]
       (map (fn [vs] (histogram-internal [vs nmn nmx intervals])) nvs))
     (histogram-internal (process-vs-and-bins vs bins-or-estimate-method mn mx)))))

;; distances

(defn- quantize-distribution
  "Find probabilities from distribution for given intervals from histogram."
  [xs distr bins]
  (let [{:keys [^double step bins]} (if (map? xs) xs (histogram xs bins))
        last-idx (dec (count bins))
        counts (map second bins)]
    [counts (map-indexed (fn [^long id [^double s]]
                           (condp = id
                             0 (r/cdf distr (+ s step))
                             last-idx (- 1.0 (r/cdf distr s))
                             (r/cdf distr s (+ s step)))) bins)]))

(defn- pq-from-histograms
  [P Q bins]
  (let [[Ph Qh] (histogram [P Q] bins)]
    [(map second (:bins Ph))
     (map second (:bins Qh))]))

(defn- normalize-PQ
  [P-observed Q-expected bins probabilities?]
  (let [[P Q] (cond (r/distribution? Q-expected) (quantize-distribution P-observed Q-expected bins)
                    bins (pq-from-histograms P-observed Q-expected bins) 
                    :else [P-observed Q-expected])]
    (cond
      (r/distribution? Q-expected) [(v/div P (v/sum P)) Q]
      probabilities? [(v/div P (v/sum P)) (v/div Q (v/sum Q))]
      :else [P Q])))

(defn- remove-zeros-pairwise
  [[P Q]]
  (let [pairs (->> (map v/vec2 P Q)
                   (remove (fn [^Vec2 v] (or (zero? (.x v))
                                            (zero? (.y v))))))]
    [(map first pairs) (map second pairs)]))

(defn- safe-div
  ^double [^double n ^double d ^double e]
  (if (zero? d) (/ n e) (/ n d)))

(defn- make-safe-log
  [^double ex ^double e]
  (cond
    (== ex m/E) (fn ^double [^double v] (m/log (if (zero? v) e v)))
    (== ex 2.0) (fn ^double [^double v] (m/log2 (if (zero? v) e v)))
    (== ex 10.0) (fn ^double [^double v] (m/log10 (if (zero? v) e v)))
    :else (fn ^double [^double v] (m/logb ex (if (zero? v) e v)))))

(defn dissimilarity
  "Various PDF distance between two histograms (frequencies) or probabilities.

  Q can be a distribution object. Then, histogram will be created out of P.

  Arguments:

  * `method` - distance method
  * `P-observed` - frequencies, probabilities or actual data (when Q is a distribution of `:bins` is set)
  * `Q-expected` - frequencies, probabilities or distribution object (when P is a data or `:bins` is set)

  Options:

  * `:probabilities?` - should P/Q be converted to a probabilities, default: `true`.
  * `:epsilon` - small number which replaces `0.0` when division or logarithm is used`
  * `:log-base` - base for logarithms, default: `e`
  * `:power` - exponent for `:minkowski` distance, default: `2.0`
  * `:bins` - number of bins or bins estimation method, see [[histogram]].

  The list of methods: `:euclidean`, `:city-block`, `:manhattan`, `:chebyshev`, `:minkowski`, `:sorensen`, `:gower`, `:soergel`, `:kulczynski`, `:canberra`, `:lorentzian`, `:non-intersection`, `:wave-hedges`, `:czekanowski`, `:motyka`, `:tanimoto`, `:jaccard`, `:dice`, `:bhattacharyya`, `:hellinger`, `:matusita`, `:squared-chord`, `:euclidean-sq`, `:squared-euclidean`, `:pearson-chisq`, `:chisq`, `:neyman-chisq`, `:squared-chisq`, `:symmetric-chisq`, `:divergence`, `:clark`, `:additive-symmetric-chisq`, `:kullback-leibler`, `:jeffreys`, `:k-divergence`, `:topsoe`, `:jensen-shannon`, `:jensen-difference`, `:taneja`, `:kumar-johnson`, `:avg`

  See more: Comprehensive Survey on Distance/Similarity Measures between Probability Density Functions by Sung-Hyuk Cha"
  (^double [method P-observed Q-expected] (dissimilarity method P-observed Q-expected nil))
  (^double [method P-observed Q-expected {:keys [bins probabilities? ^double epsilon ^double log-base ^double power remove-zeros?]
                                          :or {probabilities? true epsilon 1.0e-6 log-base m/E power 2.0}}]
   (let [pq (normalize-PQ P-observed Q-expected bins probabilities?)
         [P Q] (if remove-zeros? (remove-zeros-pairwise pq) pq)
         log (make-safe-log log-base epsilon)]
     (case method
       :euclidean (L2 P Q)
       :city-block (L1 P Q)
       :manhattan (L1 P Q)
       :chebyshev (LInf P Q)
       :minkowski (m/pow (v/sum (map (fn [^double p ^double q] (m/pow (m/abs (- p q)) power)) P Q))
                         (/ power))
       :sorensen (safe-div (L1 P Q) (+ (v/sum P) (v/sum Q)) epsilon)
       :gower (safe-div (L1 P Q) (count P) epsilon)
       :soergel (safe-div (L1 P Q) (v/sum (v/emx P Q)) epsilon)
       :kulczynski (safe-div (L1 P Q) (v/sum (v/emn P Q)) epsilon)
       :canberra (v/sum (map (fn [^double p ^double q]
                               (safe-div (m/abs (- p q)) (+ p q) epsilon)) P Q))
       :lorentzian (v/sum (map (fn [^double p ^double q] (log (m/inc (m/abs (- p q))))) P Q))
       :non-intersection (* 0.5 (L1 P Q))
       :wave-hedges (v/sum (map (fn [^double p ^double q]
                                  (safe-div (m/abs (- p q)) (m/max p q) epsilon)) P Q))
       :czekanowski (safe-div (L1 P Q) (+ (v/sum P) (v/sum Q)) epsilon)
       :motyka (safe-div (v/sum (v/emx P Q)) (+ (v/sum P) (v/sum Q)) epsilon)
       :tanimoto (let [mx (v/emx P Q)]
                   (safe-div (v/sum (v/sub mx (v/emn P Q))) (v/sum mx) epsilon))
       :jaccard (safe-div (L2sq P Q) (- (+ (v/sum (v/sq P))
                                           (v/sum (v/sq Q)))
                                        (v/sum (v/emult P Q))) epsilon)
       :dice (safe-div (L2sq P Q) (+ (v/sum (v/sq P))
                                     (v/sum (v/sq Q))) epsilon)
       :bhattacharyya (- ^double (log (v/sum (v/sqrt (v/emult P Q)))))
       :hellinger (* 2.0 (m/sqrt (- 1.0 (v/sum (v/sqrt (v/emult P Q))))))
       :matusita (m/sqrt (- 2.0 (* 2.0 (v/sum (v/sqrt (v/emult P Q))))))
       :squared-chord (L2sq (v/sqrt P) (v/sqrt Q))
       :euclidean-sq (L2sq P Q)
       :squared-euclidean (L2sq P Q)
       :pearson-chisq (v/sum (map (fn [^double p ^double q]
                                    (safe-div (m/sq (- p q)) q epsilon)) P Q))
       :chisq (v/sum (map (fn [^double p ^double q]
                            (safe-div (m/sq (- p q)) q epsilon)) P Q))
       :neyman-chisq (v/sum (map (fn [^double p ^double q]
                                   (safe-div (m/sq (- p q)) p epsilon)) P Q))
       :squared-chisq (v/sum (map (fn [^double p ^double q]
                                    (safe-div (m/sq (- p q)) (+ p q) epsilon)) P Q))
       :symmetric-chisq (* 2.0 (v/sum (map (fn [^double p ^double q]
                                             (safe-div (m/sq (- p q)) (+ p q) epsilon)) P Q)))
       :divergence (* 2.0 (v/sum (map (fn [^double p ^double q]
                                        (safe-div (m/sq (- p q)) (m/sq (+ p q)) epsilon)) P Q)))
       :clark (m/sqrt (v/sum (map (fn [^double p ^double q]
                                    (m/sq (safe-div (m/abs (- p q)) (+ p q) epsilon))) P Q)))
       :additive-symmetric-chisq (v/sum (map (fn [^double p ^double q]
                                               (safe-div (* (m/sq (- p q)) (+ p q))
                                                         (* p q)
                                                         epsilon)) P Q))
       :kullback-leibler (v/sum (map (fn [^double p ^double q]
                                       (* p ^double (log (safe-div p q epsilon)))) P Q))
       :jeffreys (v/sum (map (fn [^double p ^double q]
                               (* (- p q) ^double (log (safe-div p q epsilon)))) P Q))
       :k-divergence (v/sum (map (fn [^double p ^double q]
                                   (* p ^double (log (safe-div (* 2.0 p) (+ p q) epsilon)))) P Q))
       :topsoe (v/sum (map (fn [^double p ^double q]
                             (+ (* p ^double (log (safe-div (* 2.0 p) (+ p q) epsilon)))
                                (* q ^double (log (safe-div (* 2.0 q) (+ p q) epsilon))))) P Q))
       :jensen-shannon (* 0.5 (+ (v/sum (map (fn [^double p ^double q]
                                               (* p ^double (log (safe-div (* 2.0 p) (+ p q) epsilon)))) P Q))
                                 (v/sum (map (fn [^double p ^double q]
                                               (* q ^double (log (safe-div (* 2.0 q) (+ p q) epsilon)))) P Q))))
       :jensen-difference (v/sum (map (fn [^double p ^double q]
                                        (let [pq2 (* 0.5 (+ p q))]
                                          (- (* 0.5 (+ (* p ^double (log p))
                                                       (* q ^double (log q)))) (* pq2 ^double (log pq2))))) P Q))
       :taneja (v/sum (map (fn [^double p ^double q]
                             (let [pq2 (* 0.5 (+ p q))]
                               (* pq2 ^double (log (safe-div pq2 (m/sqrt (* p q)) epsilon))))) P Q))
       :kumar-johnson (v/sum (map (fn [^double p ^double q]
                                    (safe-div (m/sq (- (m/sq p) (m/sq q)))
                                              (* 2.0 (m/sqrt (m/cb (* p q)))) epsilon)) P Q))
       :avg (* 0.5 (+ (L1 P Q) (LInf P Q)))))))

(defn similarity
  "Various PDF similarities between two histograms (frequencies) or probabilities.

  Q can be a distribution object. Then, histogram will be created out of P.

  Arguments:

  * `method` - distance method
  * `P-observed` - frequencies, probabilities or actual data (when Q is a distribution)
  * `Q-expected` - frequencies, probabilities or distribution object (when P is a data)
  
  Options:

  * `:probabilities?` - should P/Q be converted to a probabilities, default: `true`.
  * `:epsilon` - small number which replaces `0.0` when division or logarithm is used`
  * `:bins` - number of bins or bins estimation method, see [[histogram]].

  The list of methods: `:intersection`, `:czekanowski`, `:motyka`, `:kulczynski`, `:ruzicka`, `:inner-product`, `:harmonic-mean`, `:cosine`, `:jaccard`, `:dice`, `:fidelity`, `:squared-chord`

  See more: Comprehensive Survey on Distance/Similarity Measures between Probability Density Functions by Sung-Hyuk Cha"
  (^double [method P-observed Q-expected] (similarity method P-observed Q-expected nil))
  (^double [method P-observed Q-expected {:keys [bins probabilities? ^double epsilon]
                                          :or {probabilities? true epsilon 1.0e-6}}]
   (let [[P Q] (normalize-PQ P-observed Q-expected bins probabilities?)]
     (case method
       :intersection (v/sum (v/emn P Q))
       :czekanowski (safe-div (* 2.0 (v/sum (v/emn P Q))) (+ (v/sum P) (v/sum Q)) epsilon)
       :motyka (safe-div (v/sum (v/emn P Q)) (+ (v/sum P) (v/sum Q)) epsilon)
       :kulczynski (safe-div (v/sum (v/emn P Q)) (L1 P Q) epsilon)
       :ruzicka (safe-div (v/sum (v/emn P Q)) (v/sum (v/emx P Q)) epsilon)
       :inner-product (v/dot P Q)
       :harmonic-mean (* 2.0 (v/sum (map (fn [^double p ^double q]
                                           (safe-div (* p q) (+ p q) epsilon)) P Q)))
       :cosine (safe-div (v/sum (v/emult P Q)) (* (v/mag P) (v/mag Q)) epsilon)
       ;; as in Kumar/Hassebrook paper
       #_#_:kumar-hassebrook (/ (safe-div (v/sq (v/sum (v/emult P Q)))
                                          (v/sum (v/emult (v/sq P) (v/sq Q)))
                                          epsilon)
                                (count P))
       :jaccard (let [pq (v/sum (v/emult P Q))]
                  (safe-div pq (- (+ (v/magsq P) (v/magsq Q)) pq) epsilon))
       :dice (let [pq (v/sum (v/emult P Q))]
               (safe-div (* 2.0 pq) (+ (v/magsq P) (v/magsq Q)) epsilon))
       :fidelity (v/sum (v/sqrt (v/emult P Q)))
       :squared-chord (dec (* 2.0 (v/sum (v/sqrt (v/emult P Q)))))))))

;;

(defn- weighted-variance-average
  ^double [groups]
  (reduce (fn [^double sum xs]
            (+ sum (* (dec (count xs))
                      (variance xs)))) 0.0 groups))

(defn pooled-variance
  "Calculate pooled variance for samples and method.

  Methods:
  * `:unbiased` - weighted average of variances (default)
  * `:biased` - biased version of `:unbiased`, no count correction.
  * `:avg` - average of variances"
  (^double [groups] (pooled-variance groups :unbiased))
  (^double [groups method]
   (let [agroups (map m/seq->double-array groups)]
     (case method
       :biased (/ (weighted-variance-average agroups)
                  (sum (map alength agroups)))
       :avg (/ (sum (map variance agroups))
               (count groups))
       (/ (weighted-variance-average agroups)
          (- (sum (map alength agroups)) (count groups)))))))

(defn pooled-stddev
  "Calculate pooled standard deviation for samples and method

  Methods:
  * `:unbiased` - sqrt of weighted average of variances (default)
  * `:biased` - biased version of `:unbiased`, no count correction.
  * `:avg` - sqrt of average of variances"
  (^double [groups] (m/sqrt (pooled-variance groups)))
  (^double [groups method] (m/sqrt (pooled-variance groups method))))

(defn pooled-mad
  "Calculate pooled median absolute deviation for samples.

  k is a scaling constant which equals around 1.4826 by default."
  (^double [groups] (pooled-mad groups 1.4826022185056023))
  (^double [groups ^double const]
   (let [Y (mapcat (fn [g] (let [md (median g)]
                            (v/shift g (m/- md)))) groups)]
     (m/* const (median-absolute-deviation Y)))))

;; effect size

(defn cohens-d
  "Calculate Cohen's d effect size between two groups.

  Cohen's d is a statistical measure used to quantify the magnitude of difference between two groups.

  Params:
  - group1: first sample
  - group2: second sample
  - method (optional): method for pooled stddev calculation (`:unbiased`, `:biased`, or `:avg`), defaults to `:unbiased`, see [[pooled-stddev]] for more info about methods.

  Returns effect size as double."
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
  (^double [[group1 group2]] (hedges-g group1 group2))
  (^double [group1 group2]
   (cohens-d group1 group2 :unbiased)))

(defn hedges-g-corrected
  "Cohen's d corrected for small group size"
  (^double [[group1 group2]] (hedges-g-corrected group1 group2))
  (^double [group1 group2]
   (cohens-d-corrected group1 group2 :unbiased)))

(defn hedges-g*
  "Less biased Hedges's g effect size for two groups, J term correction."
  (^double [[group1 group2]] (hedges-g* group1 group2))
  (^double [group1 group2]
   (let [df (+ (count group1) (count group2) -2)
         df2 (* 0.5 df)
         j (m/exp (- (special/log-gamma df2)
                     (m/log (m/sqrt df2))
                     (special/log-gamma (* 0.5 (dec df)))))]
     (* j (hedges-g group1 group2)))))

(defn glass-delta
  "Glass's delta effect size for two groups"
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
  (^double [[group1 group2]] (cliffs-delta group1 group2))
  (^double [group1 group2]
   (/ (sum (for [a group1
                 b group2]
             (m/signum (compare a b))))
      (* (count group1) (count group2)))))

;;

(defn ameasure
  "Vargha-Delaney A measure for two populations a and b"
  (^double [[group1 group2]] (ameasure group1 group2))
  (^double [group1 group2]
   (let [m (count group1)
         n (count group2)
         r1 (sum (take m (m/rank1 (concat group1 group2))))]
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
   (let [{kde1 :kde ^double h1 :h
          ^double mn1 :mn ^double mx1 :mx} (kd/kernel-density+ kde group1 {:bandwidth bandwidth})
         {kde2 :kde ^double h2 :h
          ^double mn2 :mn ^double mx2 :mx} (kd/kernel-density+ kde group2 {:bandwidth bandwidth})
         h (* 2.0 (+ h1 h2))
         mn (- (m/min mn1 mn2) h)
         mx (+ (m/max mx1 mx2) h)
         iters (m/max 2 min-iterations)
         ranges (partition 2 1 (m/slice-range mn mx steps))
         i1 (integrate-kde iters kde1 ranges)
         i2 (integrate-kde iters kde2 ranges)]
     (sum (map min i1 i2)))))

(defn cohens-u1-normal
  "Returns Cohen's U1 for normal samples with equal variability.

  Proportion of nonoverlap between two distributions. Symmetric.

  `method` - pooled standard deviation method: `:unbiased` (default), `:biased`, `:avg`

  Based on orignal definition by Cohen."
  (^double [group1 group2] (cohens-u1-normal group1 group2 :unbiased))
  (^double [group1 group2 method] (cohens-u1-normal (cohens-d group1 group2 method)))
  (^double [^double d]
   (let [p (r/cdf r/default-normal (m/* 0.5 (m/abs d)))]
     (m// (m/dec (m/* 2.0 p)) p))))

(defn cohens-u2-normal
  "Returns Cohen's U2 for normal samples with equal variability.

  The percentage of one group that exceeds the same percentage in the second group. Symmetric.
  
  `method` - pooled standard deviation method: `:unbiased` (default), `:biased`, `:avg`

  Based on orignal definition by Cohen."
  (^double [group1 group2] (cohens-u2-normal group1 group2 :unbiased))
  (^double [group1 group2 method] (cohens-u2-normal (cohens-d group1 group2 method)))
  (^double [^double d] (r/cdf r/default-normal (m/* 0.5 (m/abs d)))))

(defn cohens-u3-normal
  "Returns Cohen's U3 for normal samples with equal variability.

  Percentage of one group that is lower than median of the other group. Not symmetric.

  `method` - pooled standard deviation method: `:unbiased` (default), `:biased`, `:avg`

  Based on orignal definition by Cohen."
  (^double [group1 group2] (cohens-u3-normal group1 group2 :unbiased))
  (^double [group1 group2 method] (cohens-u3-normal (cohens-d group1 group2 method)))
  (^double [^double d] (r/cdf r/default-normal d)))

(defn cohens-u2
  "Returns Cohen's U2 for any samples.

  The percentage of one group that exceeds the same percentage in the second group. Symmetric."
  (^double [[group1 group2]] (cohens-u2 group1 group2))
  (^double [group1 group2]
   (let [g1 (r/distribution :real-discrete-distribution {:data group1})
         g2 (r/distribution :real-discrete-distribution {:data group2})
         target-fn (fn [^double p]
                     (let [p (m/constrain p 4.9E-324 0.9999999999999999)
                           p- (m/- 1.0 p)
                           q1 (Vec2. (r/icdf g1 p) (r/icdf g1 p-))
                           q2 (Vec2. (r/icdf g2 p-) (r/icdf g2 p))]
                       (-> (v/sub q1 q2) v/abs v/mn)))]
     (-> (opt/minimize :brent target-fn {:bounds [[0.5 1.0]]}) ffirst double))))

(defn cohens-u1
  "Returns Cohen's U1 for any samples.

  Proportion of nonoverlap between two distributions. Symmetric."
  (^double [[group1 group2]] (cohens-u1 group1 group2))
  (^double [group1 group2]
   (let [u2 (cohens-u2 group1 group2)]
     (m// (m/dec (m/* 2.0 u2)) u2))))

(defn cohens-u3
  "Returns Cohen's U3 for any samples.

  Percentage of one group that is lower than median of the other group. Not symmetric."
  (^double [[group1 group2]] (cohens-u3 group1 group2))
  (^double [group1 group2] (cohens-u3 group1 group2 :legacy))
  (^double [group1 group2 estimation-strategy]
   (let [m (median group1 estimation-strategy)]
     (/ (count (filter (fn [^double v] (m/< v m)) group2)) (double (count group2))))))

;;

(defn pearson-r
  "Pearson `r` correlation coefficient"
  (^double [[group1 group2]] (pearson-r group1 group2))
  (^double [group1 group2]
   (pearson-correlation group1 group2)))

(defn r2-determination
  "Coefficient of determination"
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
  (^double [[group1 group2]] (eta-sq group1 group2))
  (^double [group1 group2]
   (.getRSquare (local-linear-regression group1 group2))))

(defn omega-sq
  "Adjusted R2"
  (^double [[group1 group2]] (omega-sq group1 group2))
  (^double [group1 group2]
   (let [lm (local-linear-regression group1 group2)
         mse (.getMeanSquareError lm)]
     (/ (- (.getRegressionSumSquares lm) mse)
        (+ (.getTotalSumSquares lm) mse))))
  (^double [group1 group2 ^double degrees-of-freedom]
   (let [lm (local-linear-regression group1 group2)
         mse (.getMeanSquareError lm)]
     (/ (- (.getRegressionSumSquares lm) (* degrees-of-freedom mse))
        (+ (.getTotalSumSquares lm) mse)))))

(defn epsilon-sq
  "Less biased R2"
  (^double [[group1 group2]] (epsilon-sq group1 group2))
  (^double [group1 group2]
   (let [lm (local-linear-regression group1 group2)]
     (/ (- (.getRegressionSumSquares lm) (.getMeanSquareError lm))
        (.getTotalSumSquares lm)))))

(defn cohens-f2
  "Cohens f2, by default based on `eta-sq`.

  Possible `type` values are: `:eta` (default), `:omega` and `:epsilon`."
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
              [k (sum (map second v))]) m)
       (sort-by first)))

(defn contingency-table->marginals
  [ct]
  (let [rows (ct-marginals-sum (group-by ffirst ct))
        cols (ct-marginals-sum (group-by (comp second first) ct))
        n (sum (vals ct))
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
  (^double [group1 group2] (cramers-c (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (/ chi2 (+ n chi2))))))

(defn cramers-v
  "Cramer's V effect size for discrete data."
  (^double [group1 group2] (cramers-v (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long k ^long r ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (/ (/ chi2 n)
                (min (dec k) (dec r)))))))

(defn cramers-v-corrected
  "Corrected Cramer's V"
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
  (^double [group1 group2] (cohens-w (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (/ chi2 n)))))

(defn tschuprows-t
  "Tschuprows T effect size for discrete data"
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

(defn durbin-watson
  "Lag-1 Autocorrelation test for residuals"
  [rs]
  (let [es (map (fn [[^double x1 ^double x2]] (- x1 x2)) (partition 2 1 rs))]
    (/ (v/dot es es) (v/dot rs rs))))

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
    (if (every? number? xs) (map m/one? xs) xs)
    (let [f (if (sequential? true-value) (set true-value) true-value)]
      (map f xs))))

(defn- infer-confusion-matrix
  [confusion-matrix]
  (cond

    (and (map? confusion-matrix)
         (every? #{[:t :p] [:t :n] [:f :p] [:f :n]} (keys confusion-matrix)))
    (into {} (map (fn [[[a b] v]] [(keyword (str (name a) (name b))) v]) confusion-matrix))

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
  (reduce + 0.0 (map * xs1 xs2)))

;; http://feldman.faculty.pstat.ucsb.edu/174-03/lectures/l12
(defn acf
  "Calculate acf (autocorrelation function) for given number of lags or a list of lags.

  If lags is omitted function returns maximum possible number of lags.

  See also [[acf-ci]], [[pacf]], [[pacf-ci]]"
  ([data] (acf data (dec (count data))))
  ([data lags]
   (let [vdata (vec (demean data))
         rcov0 (/ (cov-for-acf vdata vdata))]
     (map (fn [^long lag]
            (if (zero? lag)
              1.0
              (let [v2 (subvec vdata lag)
                    v1 (subvec vdata 0 (count v2))]
                (* rcov0 (cov-for-acf v1 v2))))) (if (number? lags)
                                                   (range (inc (int lags)))
                                                   (seq lags))))))

;; http://feldman.faculty.pstat.ucsb.edu/174-03/lectures/l13
(defn pacf
  "Caluclate pacf (partial autocorrelation function) for given number of lags.

  If lags is omitted function returns maximum possible number of lags.

  `pacf` returns also lag `0` (which is `0.0`).

  See also [[acf]], [[acf-ci]], [[pacf-ci]]"
  ([data] (pacf data (dec (count data))))
  ([data ^long lags]
   (let [acfs (vec (acf data lags))
         phis (reductions (fn [curr ^long id]
                            (let [phi (/ (- ^double (acfs id)
                                            (sum
                                             (map-indexed (fn [^long idx ^double c]
                                                            (* c ^double (acfs (dec (- id idx))))) curr)))
                                         (- 1.0
                                            ^double (sum (map-indexed (fn [^long id ^double c]
                                                                        (* c ^double (acfs (inc id)))) curr))))]

                              (conj (mapv (fn [^double p1 ^double p2]
                                            (- p1 (* phi p2))) curr (reverse curr)) phi))) [(acfs 1)] (range 2 (inc lags)))]
     (conj (map last phis) 0.0))))

(defn- p-acf-ci-value
  ^double [data ^double alpha]
  (* (/ (m/sqrt (count data)))
     ^double (r/icdf r/default-normal (* 0.5 (inc (- 1.0 alpha))))))

(defn pacf-ci
  "[[pacf]] with added confidence interval data."
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
  ([vs] (percentile-bc-extent vs 2.5))
  ([vs ^double p] (percentile-bc-extent vs p (- 100.0 p)))
  ([vs p1 p2] (percentile-bc-extent vs p1 p2 :legacy))
  ([vs p1 p2 estimation-strategy]
   (percentile-bca-extent vs p1 p2 0.0 estimation-strategy)))

;;

(def binomial-ci-methods (sort [:asymptotic :agresti-coull :clopper-pearson :wilson :prop.test
                              :cloglog :logit :probit :arcsine]))

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

  Default alpha is 0.05
  
  Returns a triple [lower ci, upper ci, p=successes/trials]"
  ([^long number-of-successes ^long number-of-trials]
   (binomial-ci number-of-successes number-of-trials :asymptotic))
  ([^long number-of-successes ^long number-of-trials method]
   (binomial-ci number-of-successes number-of-trials method 0.05))
  ([^long number-of-successes ^long number-of-trials method ^double alpha]
   (let [p (/ (double number-of-successes) number-of-trials)
         alpha2 (* 0.5 alpha)
         ^double z (r/icdf r/default-normal (- 1.0 alpha2))
         z2 (* z z)
         x0? (zero? number-of-successes)
         xn? (== number-of-trials number-of-successes)]
     (case method
       :all (into {} (map #(vector %1 (binomial-ci number-of-successes number-of-trials % alpha))
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

(defn skewness-test
  "Normality test for skewness."
  ([xs] (skewness-test xs nil))
  ([xs params] (skewness-test xs nil params))
  ([xs skew {:keys [sides type]
             :or {sides :two-sided type :g1}}]
   (let [skew (double (or skew (skewness xs type)))
         n (count xs)
         y (* skew (m/sqrt (/ (* (inc n) (+ n 3))
                              (* 6.0 (- n 2)))))
         beta2- (dec (/ (* 3.0 (+ (* n n) (* 27 n) -70) (+ n 1) (+ n 3))
                        (* (- n 2) (+ n 5) (+ n 7) (+ n 9))))
         w2 (dec (m/sqrt (* 2.0 beta2-)))
         delta (/ 1.0 (m/sqrt (* 0.5 (m/log w2))))
         alpha (m/sqrt (/ 2.0 (dec w2)))
         ya (double (if (zero? y) (/ 1.0 alpha) (/ y alpha)))
         Z (* delta (m/log (+ ya (m/sqrt (inc (* ya ya))))))]
     {:p-value (p-value r/default-normal Z sides)
      :Z Z
      :skewness skew})))

(defn kurtosis-test
  "Normality test for kurtosis"
  ([xs] (kurtosis-test xs nil))
  ([xs params] (kurtosis-test xs nil params))
  ([xs kurt {:keys [sides type]
             :or {sides :two-sided type :kurt}}]
   (let [kurt (double (or kurt (kurtosis xs type)))
         n (count xs)
         e (/ (* 3.0 (dec n)) (inc n))
         varb2 (/ (* 24.0 (* n (- n 2) (- n 3)))
                  (* (m/sq (inc n)) (+ n 3) (+ n 5)))
         x (/ (- kurt e) (m/sqrt varb2))
         sqrtbeta1 (* (/ (* 6.0 (+ (* n n) (* -5 n) 2))
                         (* (+ n 7) (+ n 9)))
                      (m/sqrt (/ (* 6.0 (* (+ n 3) (+ n 5)))
                                 (* n (- n 2) (- n 3)))))
         a (+ 6.0 (* (/ 8.0 sqrtbeta1) (+ (/ 2.0 sqrtbeta1)
                                          (m/sqrt (inc (/ 4.0 (* sqrtbeta1 sqrtbeta1)))))))
         term1 (- 1.0 (/ 2.0 (* 9.0 a)))
         denom (inc (* x (m/sqrt (/ 2.0 (- a 4.0)))))
         term2 (* (m/signum denom) (m/cbrt (/ (- 1.0 (/ 2.0 a))
                                              (m/abs denom))))
         Z (/ (- term1 term2)
              (m/sqrt (/ 2.0 (* 9.0 a))))]
     {:p-value (p-value r/default-normal Z sides)
      :Z Z
      :kurtosis kurt})))

(defn normality-test
  "Normality test based on skewness and kurtosis"
  ([xs] (normality-test xs nil))
  ([xs params] (normality-test xs nil nil params))
  ([xs skew kurt {:keys [sides]
                  :or {sides :one-sided-greater}}]
   (let [{^double skew-Z :Z skew :skewness} (skewness-test xs skew nil)
         {^double kurt-Z :Z kurt :kurtosis} (kurtosis-test xs kurt nil)
         Z (+ (* skew-Z skew-Z)
              (* kurt-Z kurt-Z))]
     {:p-value (p-value (r/distribution :chi-squared {:degrees-of-freedom 2}) Z sides)
      :Z Z
      :skewness skew
      :kurtosis kurt})))

(defn jarque-bera-test
  "Goodness of fit test whether skewness and kurtosis of data match normal distribution"
  ([xs] (jarque-bera-test xs nil))
  ([xs params] (jarque-bera-test xs nil nil params))
  ([xs skew kurt {:keys [sides]
                  :or {sides :one-sided-greater}}]
   (let [skew (double (or skew (skewness xs :g1)))
         kurt (double (or kurt (kurtosis xs :g2)))
         n (count xs)
         Z (* n m/SIXTH (+ (* skew skew) (* 0.25 kurt kurt)))]
     {:p-value (p-value (r/distribution :chi-squared {:degrees-of-freedom 2}) Z sides)
      :Z Z
      :skewness skew
      :kurtosis (+ kurt 3.0)})))

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
      :confidence-interval (let [bci (partial binomial-ci number-of-successes number-of-trials ci-method)]
                             (sides-case sides
                                         (vec (butlast (bci (- 1.0 alpha))))
                                         [(first (bci (- 1.0 (* alpha 2.0)))) 1.0]
                                         [0.0 (second (bci (- 1.0 (* alpha 2.0))))]))})))

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
       (-> (t-test-one-sample (map - xs ys) params)
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
       (-> (z-test-one-sample (map - xs ys) params)
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
  "Goodness of fit"
  [xs p ^double lambda]
  (let [cnt (count xs)
        n (sum xs)
        df (dec cnt)
        p (or p (repeat cnt 1.0))
        psum (sum p)
        p (map (fn [^double p] (/ p psum)) p)
        xhat (map (fn [^double p] (* n p)) p)
        stat (condp = (double lambda)
               0.0 (* 2.0 (sum (map (fn [^long a ^double b]
                                      (* a (- (m/log a) (m/log b)))) xs xhat)))
               -1.0 (* 2.0 (sum (map (fn [^double a ^long b]
                                       (* a (- (m/log a) (m/log b)))) xhat xs)))
               (* (/ 2.0 (* lambda (inc lambda)))
                  (sum (map (fn [^long a ^double b]
                              (* a (dec (m/pow (/ a b) lambda)))) xs xhat))))]
    {:stat stat :df df :n n :expected xhat :p p
     :estimate (map (fn [^double v] (/ v n)) xs)}))

(defn- pdt-distribution
  [xs distr bins ^double lambda]
  (let [[counts probabilities] (quantize-distribution xs distr bins)]
    (pdt-gof counts probabilities lambda)))

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
  "Power divergence test.

  First argument should be one of:
  
  * contingency table
  * sequence of counts (for goodness of fit)
  * sequence of data (for goodness of fit against distribution)

  For goodness of fit there are two options:
  
  * comparison of observed counts vs expected probabilities or weights (`:p`)
  * comparison of data against given distribution (`:p`), in this case histogram from data is created and compared to distribution PDF in bins ranges. Use `:bins` option to control histogram creation.

  Options are:
  
  * `:lambda` - test type:
      * `1.0` - [[chisq-test]]
      * `0.0` - [[multinomial-likelihood-ratio-test]]
      * `-1.0` - [[minimum-discrimination-information-test]]
      * `-2.0` - [[neyman-modified-chisq-test]]
      * `-0.5` - [[freeman-tukey-test]]
      * `2/3` - [[cressie-read-test]] - default
  * `:p` - probabilites, weights or distribution object.
  * `:alpha` - significance level (default: 0.05)
  * `:ci-sides` - confidence interval sides (default: `:two-sided`)
  * `:sides`  - p-value sides (`:two-sided`, `:one-side-greater` - default, `:one-side-less`)
  * `:bootstrap-samples` - number of samples to estimate confidence intervals (default: 1000)
  * `:ddof` - delta degrees of freedom, adjustment for dof (default: 0.0)
  * `:bins` - number of bins or estimator name for histogram"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {}))
  ([contingency-table-or-xs {:keys [^double lambda ci-sides sides p ^double alpha ^long bootstrap-samples
                                    ^long ddof bins]
                             :or {lambda m/TWO_THIRD sides :one-sided-greater ci-sides :two-sided
                                  alpha 0.05 bootstrap-samples 1000 ddof 0}}]
   (let [{:keys [df stat] :as res} (-> (cond
                                         (and p (r/distribution? p))
                                         (pdt-distribution contingency-table-or-xs p bins lambda)

                                         (and (sequential? contingency-table-or-xs)
                                              (every? number? contingency-table-or-xs))
                                         (pdt-gof contingency-table-or-xs p lambda)

                                         :else (pdt-multi contingency-table-or-xs lambda))
                                       (update :df (fn [^long df] (- df ddof))))
         distr (r/distribution :chi-squared {:degrees-of-freedom df})
         res (assoc res :lambda lambda :sides sides :test-type sides :ci-sides ci-sides :chi2 stat :alpha alpha :level (- 1.0 alpha)
                    :p-value (p-value distr stat sides))]
     (assoc res :confidence-interval (pdt-bootstrap-ci res bootstrap-samples)))))

(defn chisq-test
  "Chi square test, a power divergence test for `lambda` 1.0"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda 1.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda 1.0))))

(defn multinomial-likelihood-ratio-test
  "Multinomial likelihood ratio test, a power divergence test for `lambda` 0.0"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda 0.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda 0.0))))

(defn minimum-discrimination-information-test
  "Minimum discrimination information test, a power divergence test for `lambda` -1.0"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda -1.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda -1.0))))

(defn neyman-modified-chisq-test
  "Neyman modifield chi square test, a power divergence test for `lambda` -2.0"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda -2.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda -2.0))))

(defn freeman-tukey-test
  "Freeman-Tukey test, a power divergence test for `lambda` -0.5"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda -0.5}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda -0.5))))

(defn cressie-read-test
  "Cressie-Read test, a power divergence test for `lambda` 2/3"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda m/TWO_THIRD}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda m/TWO_THIRD))))

;; copy docs
(doseq [v [#'chisq-test #'multinomial-likelihood-ratio-test #'minimum-discrimination-information-test
           #'neyman-modified-chisq-test #'freeman-tukey-test #'cressie-read-test]]
  (alter-meta! v update :doc str "\n\n" (:doc (meta #'power-divergence-test))))

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
                        (reductions m/+ 0)
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
    (reduce - (- n) (map (fn [^long idx]
                           (* (/ (+ idx idx 1.0) n)
                              (+ (m/log (r/cdf d (aget xs idx)))
                                 (m/log (r/ccdf d (aget xs (- n idx 1))))))) (range n)))))

(defn ad-test-one-sample
  "Performs the Anderson-Darling (AD) test to assess whether a sample comes from a specified theoretical distribution.

  Parameters:
  - `xs` (seq of numbers): Sample data to be tested.
  - `distribution-or-ys`: Distribution object (default: normal) or a second dataset for an empirical comparison. In case of dataset, kernel density destination (kde) is performed. 
  - Options:
    - `:sides` (keyword, default `:right`): Specifies the alternative hypothesis.
      - `:two-sided` tests for any difference.
      - `:right` (default) tests if `xs` is stochastically greater.
      - `:left` tests if `xs` is stochastically smaller.
    - `:kernel` (keyword, default `:gaussian`): Kernel method for density estimation. When kernel is set to `:enumerated`, enumerated distribution is created instead. 
    - `:bandwidth` (double, optional): Bandwidth for kernel density estimation.

  Returns a map containing:
    - `:stat`: AD test statistic.
    - `:A2`: Anderson-Darling test statistic.
    - `:mean`: Mean of the sample.
    - `:stddev`: Standard deviation of the sample.
    - `:n`: Sample size.
    - `:sides`: Alternative hypothesis used.
    - `:p-value`: Probability of observing the result under the null hypothesis."
  ([xs] (ad-test-one-sample xs r/default-normal))
  ([xs distribution-or-ys] (ad-test-one-sample xs distribution-or-ys {}))
  ([xs distribution-or-ys {:keys [sides kernel bandwidth] :or {sides :right kernel :gaussian}}]
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
  "Performs the Kolmogorov-Smirnov (KS) test to compare a sample distribution against a theoretical distribution or another empirical sample.

  Parameters:
  - `xs` (seq of numbers): Sample data to be tested.
  - `distribution-or-ys`: Distribution object (default: normal) or a second dataset for an empirical comparison. In case of dataset, kernel density destination (kde) is performed. 
  - Options:
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis.
      - `:two-sided` (default) tests for any difference.
      - `:right` tests if `xs` is stochastically greater.
      - `:left` tests if `xs` is stochastically smaller.
    - `:kernel` (keyword, default `:gaussian`): Kernel method for density estimation. When kernel is set to `:enumerated`, enumerated distribution is created instead. 
    - `:bandwidth` (double, optional): Bandwidth for kernel density estimation.
    - `:distinct?` (boolean, default `true`): Whether to remove duplicate values before computation.

  Returns:
  - A map containing:
    - `:n`: Sample size.
    - `:dp`: Maximum positive difference between empirical and reference CDF.
    - `:dn`: Maximum negative difference.
    - `:d`: KS test statistic (max absolute difference).
    - `:stat`: Scaled KS test statistic.
    - `:p-value`: Probability of observing the result under the null hypothesis."
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
         ^double dp (reduce max (map - (rest idxs) cdfs))
         dn (- ^double (reduce min (map - (butlast idxs) cdfs)))
         d (max dp dn)]
     {:n n :dp dp :dn dn :d d :sides sides
      :stat (sides-case sides d dp dn) 
      :p-value (sides-case sides
                           (p-value (r/distribution :kolmogorov-smirnov {:n n}) d :right)
                           (p-value (r/distribution :kolmogorov-smirnov+ {:n n}) dp :right)
                           (p-value (r/distribution :kolmogorov-smirnov+ {:n n}) dn :right))})))

(defn- process-ks-diffs
  [vs os ^long nx ^long ny]
  (let [cnt- (dec (count vs))
        dx (/ 1.0 nx)
        dy (/ -1.0 ny)]
    (loop [i (long 0)
           d (double 0.0)
           dn (double 0.0)
           dp (double 0.0)]
      (let [id (long (os i))
            nd (+ d (if (< id nx) dx dy))]
        (if (= i cnt-)
          [(min nd dn) (max nd dp)]
          (let [v1 (double (vs id))
                v2 (double (vs (os (inc i))))]
            (if (m/not== v1 v2)
              (recur (inc i) nd (min nd dn) (max nd dp))
              (recur (inc i) nd dn dp))))))))

;; reimplementing ACM exact version using https://arxiv.org/pdf/2102.08037
(defn- ks-exact
  [^double d ^long m ^long n]
  (let [cnm (unchecked-long (m/ceil (m/* (m/- d 1.0e-12) n m)))
        lag (double-array (map (fn [^long k]
                                 (if (m/>= (m/* m (m/inc k)) cnm) 1.0 0.0)) (range n)))]
    (loop [k (long 1)
           lst (double 0.0)]
      (if (m/> k m)
        lst
        (recur (m/inc k) (double (loop [l (long 1)
                                        lst (if (m/>= (m/* k n) cnm) 1.0 0.0)]
                                   (if (m/> l n)
                                     lst
                                     (let [l- (m/dec l)
                                           v (if (m/>= (m/abs (m/- (m/* k n) (m/* l m))) cnm)
                                               1.0
                                               (m// (m/+ (m/* k (Array/aget lag l-))
                                                         (m/* l lst))
                                                    (m/+ k l)))]
                                       (Array/aset lag l- v)
                                       (recur (m/inc l) v))))))))))

(defn ks-test-two-samples
  "Performs the two-sample Kolmogorov-Smirnov (KS) test to compare the distributions of two independent samples, `xs` and `ys`. This test determines whether the two samples come from the same distribution.

  Arguments:
  - `xs`: First sample (a sequence of numerical values).
  - `ys`: Second sample (a sequence of numerical values).
  - options map:
    - `:method`: `:exact` or `approximate` (asymptotic, default)
    - `:sides` (default: `:two-sided`): Specifies the alternative hypothesis. 
      Possible values:
      - `:two-sided`: The two distributions are different.
      - `:right`: The distribution of `xs` is stochastically greater than `ys`.
      - `:left`: The distribution of `xs` is stochastically less than `ys`.
    - `:distinct?` (default: `true`): If true, removes duplicate values from `xs` and `ys` before computation.

  Returns:
  A map with the following keys:
  - `:n`: Effective sample size.
  - `:nx`: Number of observations in `xs`.
  - `:ny`: Number of observations in `ys`.
  - `:dp`: Maximum positive difference between empirical cumulative distribution functions.
  - `:dn`: Maximum negative difference between empirical cumulative distribution functions.
  - `:d`: The maximum absolute difference between ECDFs.
  - `:stat`: KS statistic, scaled `d` for asymptotic method.
  - `:KS`: Alias for `:stat`.
  - `:sides`: The alternative hypothesis used.
  - `:p-value`: The computed p-value indicating the significance of the test.

  Please note that `:right` or `:left` sides can give wrong results for `:exact` method."
  ([xs ys] (ks-test-two-samples xs ys {}))
  ([xs ys {:keys [method sides distinct?] :or {method :approximate sides :two-sided distinct? true}}]
   (let [xs (if distinct? (distinct xs) xs)
         ys (if distinct? (distinct ys) ys)
         nx (count xs)
         ny (count ys)
         vs (vec (concat xs ys))
         sort-idxs (vec (m/order vs))         
         [dn dp] (process-ks-diffs vs sort-idxs nx ny)
         dn (- (double dn))
         dp (double dp)
         d (max dn dp)
         res {:nx nx :ny ny :dp dp :dn dn :d d
              :method method
              :sides sides}]
     (if (= method :exact)
       (let [n (+ nx ny)
             stat (sides-case sides d dp dn)]
         (assoc res :n n :stat stat :KS stat
                :p-value (ks-exact stat nx ny)))
       (let [n (/ (* nx ny) (double (+ nx ny)))
             stat (* (m/sqrt n) (sides-case sides d dp dn))]
         (assoc res :n n :stat stat :KS stat
                :p-value (sides-case sides
                                     (p-value (r/distribution :kolmogorov) stat :right)
                                     (m/exp (* -2.0 (* stat stat)))
                                     (m/exp (* -2.0 (* stat stat))))))))))

(defn kruskal-test
  "Performs the Kruskal-Wallis H-test (rank sum test) for independent samples.

  The Kruskal-Wallis test is a non-parametric alternative to one-way ANOVA.
  It determines whether there is a statistically significant difference between the distributions of two or more independent groups. It does not assume normality but requires that distributions have a similar shape for the test to be valid.

  Parameters:

  - `data-groups` (vector of sequences): A collection where each element is a sequence 
    representing a group of observations.
  - a map containing `:sides` key with values of: `:right` (default), `:left` or `:both`

  Returns a map containing:

  - `:stat`: The Kruskal-Wallis H statistic.
  - `:n`: Total number of observations across all groups.
  - `:df`: Degrees of freedom (number of groups - 1).
  - `:k`: Number of groups.
  - `:sides`: Test side
  - `:p-value`: The p-value for the test (null hypothesis: all groups have the same distribution)."
  ([xss] (kruskal-test xss {}))
  ([xss {:keys [sides] :or {sides :right}}] ;; as in R
   (let [k (count xss)
         df (dec k)
         xs (flatten xss)
         groups (mapcat (fn [[xs id]] (repeat (count xs) id)) (map vector xss (range)))
         n (count xs)
         r (m/rank1 xs)
         ties (vals (frequencies xs))
         stat (->> (map vector r groups)
                   (group-by second)
                   (vals)
                   (map #(let [ranks (map first %)]
                           (/ (m/sq (sum ranks)) (count ranks))))
                   (sum))
         stat (/ (- (/ (* 12.0 stat)
                       (* n (inc n)))
                    (* 3.0 (inc n)))
                 (- 1.0 (/ (sum (map (fn [^long t] (- (m/cb t) t)) ties))
                           (- (m/cb n) n))))]
     {:stat stat :n n :df df :k k :sides sides
      :p-value (p-value (r/distribution :chi-squared {:degrees-of-freedom df}) stat sides)})))

;; transformations

(defn- box-cox-scaled
  ([nxs sgn ^double lambda] (box-cox-scaled nxs sgn lambda (geomean nxs)))
  ([nxs sgn ^double lambda ^double gm]
   (if (m/zero? lambda)
     (map (fn [^double s ^double x] (m/* gm s (m/log x))) sgn nxs)
     (let [fact (m/* lambda (m/pow gm (m/dec lambda)))]
       (map (fn [^double s ^double x] (m// (m/dec (m/* s (m/pow x lambda))) fact)) sgn nxs)))))

(defn- box-cox-not-scaled
  [nxs sgn ^double lambda]
  (if (m/zero? lambda)
    (map (fn [^double s ^double x] (if (m/zero? x) 0.0 (m/* s (m/log x)))) sgn nxs)
    (map (fn [^double s ^double x] (m// (m/dec (m/* s (m/pow x lambda))) lambda)) sgn nxs)))

(defn- box-cox-prepare-data
  [xs {:keys [^double alpha negative?]
       :or {alpha 0.0}}]
  [(cond-> xs
     (not (m/zero? alpha)) (v/shift alpha)
     negative? (v/abs))
   (if negative? (map m/signum xs) (repeat 1.0))])

(defn- box-cox-maximize-ll
  ^double [nxs sgn lambda-range]
  (let [[^double lambda-min ^double lambda-max :as lr] (or lambda-range [-3.0 3.0])
        lxs (v/sum (map m/log nxs))
        n- (m/- (m/* 0.5 (count nxs)))
        target (fn ^double [^double l]
                 (let [res (box-cox-not-scaled nxs sgn l)
                       v (variance res)]
                   (m/+ (m/* n- (m/log v))
                        (m/* (m/dec l) lxs))))]
    (-> (lbfgsb/maximize target {:bounds [lr]
                                 :initial [(m/lerp lambda-min lambda-max 0.51)]})
        (ffirst))))

(defn box-cox-infer-lambda
  "Find optimal `lambda` parameter for Box-Cox tranformation using maximum log likelihood method."
  (^double [xs] (box-cox-infer-lambda xs nil))
  (^double [xs lambda-range] (box-cox-infer-lambda xs lambda-range nil))
  (^double [xs lambda-range opts]
   (let [[nxs sgn] (box-cox-prepare-data xs opts)]
     (box-cox-maximize-ll nxs sgn lambda-range))))

(defn box-cox-transformation
  "Applies Box-Cox transformation to a data.

   The Box-Cox transformation is a family of power transformations used to stabilize variance and make data more normally distributed.

  Parameters:
  - `xs` (seq of numbers): The input data.
  - `lambda` (default `0.0`): The power parameter. If `nil` or `[lambda-min, lambda-max]`, `lambda` is inferred using maximum log likelihood.
  - Options map:
    - `alpha` (optional): A shift parameter applied before transformation.
    - `scaled?` (default `false`): scale by geometric mean
    - `negative?` (default `false`): allow negative values

  Returns transformed data.

  Related: `yeo-johnson-transformation`"
  ([xs] (box-cox-transformation xs nil))
  ([xs lambda] (box-cox-transformation xs lambda nil))
  ([xs lambda {:keys [scaled?] :as opts}]
   (let [[nxs sgn] (box-cox-prepare-data xs opts)
         lambda (if (number? lambda) lambda (box-cox-maximize-ll nxs sgn lambda))]
     (if scaled?
       (if (number? scaled?)
         (box-cox-scaled nxs sgn lambda scaled?)
         (box-cox-scaled nxs sgn lambda))
       (box-cox-not-scaled nxs sgn lambda)))))

(defn- yeo-johnson
  [nxs ^double lambda]
  (let [l2 (- 2.0 lambda)]
    (map (fn [^double x]
           (if (m/neg? x)
             (if (== lambda 2.0)
               (- (m/log (- 1.0 x)))
               (- (/ (dec (m/pow (- 1.0 x) l2)) l2)))
             (if (m/zero? lambda)
               (m/log (inc x))
               (/ (dec (m/pow (inc x) lambda)) lambda)))) nxs)))

(defn- yeo-johnson-maximize-ll
  ^double [nxs lambda-range]
  (let [[^double lambda-min ^double lambda-max :as lr] (or lambda-range [-3.0 3.0])
        lxs (v/sum (map (fn [^double x] (m/* (m/signum x) (m/log (m/inc (m/abs x))))) nxs))
        n- (m/- (m/* 0.5 (count nxs)))
        target (fn ^double [^double l]
                 (let [res (yeo-johnson nxs l)
                       v (variance res)]
                   (m/+ (m/* n- (m/log v))
                        (m/* (m/dec l) lxs))))]
    (-> (lbfgsb/maximize target {:bounds [lr]
                                 :initial [(m/lerp lambda-min lambda-max 0.51)]})
        (ffirst))))

(defn yeo-johnson-infer-lambda
  "Find optimal `lambda` parameter for Yeo-Johnson tranformation using maximum log likelihood method."
  (^double [xs] (yeo-johnson-infer-lambda xs nil))
  (^double [xs lambda-range] (yeo-johnson-infer-lambda xs lambda-range nil))
  (^double [xs lambda-range {:keys [^double alpha] :or {alpha 0.0}}]
   (let [nxs (if (m/zero? alpha) xs (v/shift xs alpha))]
     (yeo-johnson-maximize-ll nxs lambda-range))))

(defn yeo-johnson-transformation
  "Applies the Yeo-Johnson transformation to a dataset.

  This transformation is used to stabilize variance and make data more normally distributed. It extends the Box-Cox transformation to allow for zero and negative values.

  Parameters:
  - `xs`: The input dataset.
  - `lambda` (default: 0.0): The power parameter controlling the transformation. If `lambda` is `nil` or a range `[lambda-min, lambda-max]` it will be inferred using maximum log-likelihood method.
  - Options map:
    - `alpha` (optional): A shift parameter applied before transformation.

  Returns:
  - A transformed sequence of numbers.

  Related: `box-cox-tranformation`"
  ([xs] (yeo-johnson-transformation xs nil))
  ([xs lambda] (yeo-johnson-transformation xs lambda nil))
  ([xs lambda {:keys [^double alpha] :or {alpha 0.0}}]
   (let [nxs (if (m/zero? alpha) xs (v/shift xs alpha))
         lambda (if (number? lambda) lambda (yeo-johnson-maximize-ll nxs lambda))]
     (yeo-johnson nxs lambda))))

;;

(defn power-transformation
  "Applies a power transformation to a data."
  {:deprecated "Use `(box-cox-transformation xs lambda {:scaled true})"}
  ([xs] (power-transformation xs 0.0))
  ([xs ^double lambda]
   (box-cox-transformation xs lambda {:scaled? true}))
  ([xs ^double lambda ^double alpha]
   (box-cox-transformation xs lambda {:scaled? true :alpha alpha})))

(defn modified-power-transformation
  "Applies a modified power transformation (Bickel and Doksum) to a data."
  {:deprecated "Use `(box-cox-transformation xs lambda {:negative? true})"}
  ([xs] (modified-power-transformation xs 0.0))
  ([xs ^double lambda]
   (box-cox-transformation xs lambda {:negative? true}))
  ([xs ^double lambda ^double alpha]
   (box-cox-transformation xs lambda {:negative? true :alpha alpha})))



(m/unuse-primitive-operators)

