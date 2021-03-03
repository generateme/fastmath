(ns fastmath.stats
  "Statistics functions.

  * Descriptive statistics.
  * Correlation / covariance
  * Outliers
  * Confidence intervals
  * Extents
  * Effect size
  * Student's t-test
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
            [fastmath.random :as r])
  (:import [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.stat.descriptive.rank Percentile Percentile$EstimationType]
           [org.apache.commons.math3.stat.descriptive.moment Kurtosis Skewness]
           [org.apache.commons.math3.stat.correlation KendallsCorrelation SpearmansCorrelation PearsonsCorrelation]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn mode
  "Find the value that appears most often in a dataset `vs`.

  See also [[modes]]."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [m (StatUtils/mode (m/seq->double-array vs))]
    (aget ^doubles m 0)))

(defn modes
  "Find the values that appears most often in a dataset `vs`.

  Returns sequence with all most appearing values in increasing order.

  See also [[mode]]."
  {:metadoc/categories #{:stat}}
  [vs]
  (seq ^doubles (StatUtils/mode (m/seq->double-array vs))))

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
   (let [^Percentile perc (.withEstimationType ^Percentile (Percentile.) (or (estimation-strategies-list estimation-strategy) Percentile$EstimationType/LEGACY))]
     (.evaluate perc (m/seq->double-array vs) p))))

(defn percentiles
  "Calculate percentiles of a `vs`.

  Percentiles are sequence of values from range 0-100.

  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html).

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html)

  See also [[quantile]]."
  {:metadoc/categories #{:stat}}
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
  ([vs qs]
   (percentiles vs (map #(m/constrain (* ^double % 100.0) 0.0 100.0) qs)))
  ([vs qs estimation-strategy]
   (percentiles vs (map #(m/constrain (* ^double % 100.0) 0.0 100.0) qs) estimation-strategy)))

(defn median
  "Calculate median of `vs`. See [[median-3]]."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (percentile vs 50.0))

(defn median-3
  "Median of three values. See [[median]]."
  {:metadoc/categories #{:stat}}
  ^double [^double a ^double b ^double c]
  (m/max (m/min a b) (m/min (m/max a b) c)))

(defn mean
  "Calculate mean of `vs`"
  {:metadoc/categories #{:stat}}
  (^double [vs] (StatUtils/mean (m/seq->double-array vs))))

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

(defn median-absolute-deviation
  "Calculate MAD"
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [m (median vs)]
    (median (map (fn [^double x] (m/abs (- x m))) vs))))

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
         q1 (percentile avs 25.0 estimation-strategy)
         q3 (percentile avs 75.0 estimation-strategy)]
     (adjacent-values avs q1 q3)))
  ([vs ^double q1 ^double q3]
   (let [avs (m/seq->double-array vs)
         iqr (* 1.5 (- q3 q1))
         lav-thr (- q1 iqr)
         uav-thr (+ q3 iqr)]
     (java.util.Arrays/sort avs)
     [(first (filter #(>= (double %) lav-thr) avs))
      (last (filter #(<= (double %) uav-thr) avs))
      (median vs)])))

(defn outliers
  "Find outliers defined as values outside outer fences.

  Let Q1 is 25-percentile and Q3 is 75-percentile. IQR is `(- Q3 Q1)`.

  * LIF (Lower Outer Fence) equals `(- Q1 (* 1.5 IQR))`.
  * UIF (Upper Outer Fence) equals `(+ Q3 (* 1.5 IQR))`.

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
   (let [;;avs (m/seq->double-array vs)
         iqr (* 1.5 (- q3 q1))
         lof-thr (- q1 iqr)
         uof-thr (+ q3 iqr)]
     ;; (java.util.Arrays/sort avs)
     (filter #(let [v (double %)]
                (or (< v lof-thr)
                    (> v uof-thr))) vs))))

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

(defn kurtosis
  "Calculate kurtosis from sequence."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [^Kurtosis k (Kurtosis.)]
    (.evaluate k (m/seq->double-array vs))))

(defn moment
  "Calculate moment (central or/and absolute) of given order (default: 2).

  Additional parameters as a map:

  * `:absolute?` - calculate sum as absolute values (default: `false`)
  * `:mean?` - returns mean (proper moment) or just sum of differences (default: `true`)
  * `:center` - value of central (default: `nil` = mean)"
  (^double [vs] (moment vs 2.0 nil))
  (^double [vs ^double order] (moment vs order nil))
  (^double [vs ^double order {:keys [absolute? center mean?]
                              :or {absolute? false center nil mean? true}}]
   (let [in (double-array vs)
         cin (alength in)
         ^double center (or center (mean in))
         f (case order
             1.0 m/fast-identity
             2.0 m/sq
             3.0 m/cb
             4.0 (fn ^double [^double diff] (m/sq (m/sq diff)))
             (fn ^double [^double diff] (m/pow diff order)))
         a (if absolute? m/abs m/fast-identity)]
     (loop [idx (int 0)]
       (when (< idx cin)
         (aset in idx ^double (f (a (- (aget in idx) center))))
         (recur (inc idx))))
     (if mean? (mean in) (sum in)))))

(def ^{:deprecated "Use `moment` function"} second-moment moment)

(defn skewness
  "Calculate kurtosis from sequence."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [^Skewness k (Skewness.)]
    (.evaluate k (m/seq->double-array vs))))

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
         deltas (m/seq->double-array (repeatedly samples #(- ^double (stat-fn (r/->seq dist cnt)) m)))
         q1 (quantile deltas alpha)
         q2 (quantile deltas (- 1.0 alpha))]
     [(- m q1) (- m q2) m])))

(defn bootstrap
  "Generate set of samples of given size from provided data.

  Default `samples` is 50, number of `size` defaults to 1000"
  ([vs] (bootstrap vs 50))
  ([vs samples] (bootstrap vs samples 1000))
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
         [lav uav] (adjacent-values avs q1 q3)]
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

(stats-map [1 2 3 4 5 11])

(defn standardize
  "Normalize samples to have mean = 0 and stddev = 1."
  {:metadoc/categories #{:norm}}
  [vs]
  (seq ^doubles (StatUtils/normalize (m/seq->double-array vs))))

(defn demean
  "Subtract mean from sequence"
  {:metadoc/categories #{:norm}}
  [vs]
  (let [m (mean vs)]
    (map (fn [^double v]
           (- v m)) vs)))

(defn covariance
  "Covariance of two sequences."
  {:metadoc/categories #{:corr}}
  ^double [vs1 vs2]
  (smile.math.MathEx/cov (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn covariance-matrix
  "Generate covariance matrix from seq of seqs. Row order."
  {:metadoc/categories #{:corr}}
  [vss]
  (let [avss (map-indexed (fn [id v] [id (m/seq->double-array v)]) vss)
        cache (atom {})]
    (for [[id1 ^doubles a] avss]
      (mapv (fn [[id2 ^doubles b]]
              (let [key (sort [id1 id2])]
                (if (contains? @cache key)
                  (@cache key)
                  (let [cov (smile.math.MathEx/cov a b)]
                    (swap! cache assoc key cov)
                    cov)))) avss))))

(defn correlation
  "Correlation of two sequences."
  {:metadoc/categories #{:corr}}
  ^double [vs1 vs2]
  (smile.math.MathEx/cor (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn spearman-correlation
  "Spearman's correlation of two sequences."
  {:metadoc/categories #{:corr}}
  ^double [vs1 vs2]
  (.correlation ^SpearmansCorrelation (SpearmansCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn pearson-correlation
  "Pearson's correlation of two sequences."
  {:metadoc/categories #{:corr}}
  ^double [vs1 vs2]
  (.correlation ^PearsonsCorrelation (PearsonsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn kendall-correlation
  "Kendall's correlation of two sequences."
  {:metadoc/categories #{:corr}}
  ^double [vs1 vs2]
  (.correlation ^KendallsCorrelation (KendallsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn kullback-leibler-divergence
  "Kullback-Leibler divergence of two sequences."
  {:metadoc/categories #{:corr}}
  ^double [vs1 vs2]
  (smile.math.MathEx/KullbackLeiblerDivergence (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn jensen-shannon-divergence
  "Jensen-Shannon divergence of two sequences."
  {:metadoc/categories #{:corr}}
  ^double [vs1 vs2]
  (smile.math.MathEx/JensenShannonDivergence (m/seq->double-array vs1) (m/seq->double-array vs2)))

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

  Possible methods are: `:sqrt` `:sturges` `:rice` `:doane` `:scott` `:freedman-diaconis` (default)."
  {:metadoc/categories #{:stat}}
  ([vs] (estimate-bins vs :freedman-diaconis))
  ([vs bins-or-estimate-method]
   (if-not (keyword? bins-or-estimate-method)
     (or bins-or-estimate-method (estimate-bins vs))
     (let [n (count vs)]
       (int (condp = bins-or-estimate-method
              :sqrt (m/sqrt n)
              :sturges (inc (m/ceil (m/log2 n)))
              :rice (m/ceil (* 2.0 (m/cbrt n)))
              :doane (+ (inc (m/log2 n))
                        (m/log2 (inc (/ (m/abs (skewness vs))
                                        (m/sqrt (/ (* 6.0 (- n 2.0))
                                                   (* (inc n) (+ n 3.0))))))))
              :scott (let [vvs (m/seq->double-array vs)
                           h (/ (* 3.5 (stddev vs))
                                (m/cbrt n))]
                       (scott-fd-helper vvs h))
              (let [vvs (m/seq->double-array vs)
                    h (/ (* 2.0 (iqr vvs))
                         (m/cbrt n))]
                (scott-fd-helper vvs h))))))))

(defn histogram
  "Calculate histogram.

  Returns map with keys:

  * `:size` - number of bins
  * `:step` - distance between bins
  * `:bins` - list of pairs of range lower value and number of hits
  * `:min` - min value
  * `:max` - max value
  * `:samples` - number of used samples

  For estimation methods check [[estimate-bins]]."
  {:metadoc/categories #{:stat}}
  ([vs] (histogram vs :freedman-diaconis))
  ([vs bins-or-estimate-method] (histogram vs (estimate-bins vs bins-or-estimate-method) (extent vs)))
  ([vs ^long bins [^double mn ^double mx]]
   (let [mx (m/next-double mx)
         vs (filter #(<= mn ^double % mx) vs)
         diff (- mx mn)
         step (/ diff bins)
         search-array (double-array (map #(+ mn (* ^long % step)) (range bins)))
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
;;;;;;;;;;;;;;
;; tests


(defn- cohens-d-with-correct
  "Cohen's d effect size for two groups"
  ^double [group1 group2 ^double correct]
  (let [group1 (m/seq->double-array group1)
        group2 (m/seq->double-array group2)
        diff (- (mean group1) (mean group2))
        var1 (variance group1)
        var2 (variance group2)
        n1 (alength group1)
        n2 (alength group2)
        pooled-var (/ (+ (* n1 var1) (* n2 var2)) (+ n1 n2 correct))]
    (/ diff (m/sqrt pooled-var))))

(defn cohens-d-orig
  "Original version of Cohen's d effect size for two groups"
  {:metadoc/categories #{:effect}}
  ^double [group1 group2]
  (cohens-d-with-correct group1 group2 -2.0))

(defn cohens-d
  "Cohen's d effect size for two groups"
  {:metadoc/categories #{:effect}}
  ^double [group1 group2]
  (cohens-d-with-correct group1 group2 0.0))

(defn glass-delta
  "Glass's delta effect size for two groups"
  {:metadoc/categories #{:effect}}
  ^double [group1 group2]
  (let [group2 (m/seq->double-array group2)]
    (/ (- (mean group1) (mean group2)) (stddev group2))))

(defn hedges-g
  "Hedges's g effect size for two groups"
  {:metadoc/categories #{:effect}}
  ^double [group1 group2]
  (let [group1 (m/seq->double-array group1)
        group2 (m/seq->double-array group2)
        diff (- (mean group1) (mean group2))
        var1 (variance group1)
        var2 (variance group2)
        n1 (dec (alength group1))
        n2 (dec (alength group2))
        pooled-var (/ (+ (* n1 var1) (* n2 var2)) (+ n1 n2 -2.0))]
    (/ diff (m/sqrt pooled-var))))

(defn hedges-g*
  "Less biased Hedges's g effect size for two groups"
  {:metadoc/categories #{:effect}}
  ^double [group1 group2]
  (let [j (- 1.0 (/ 3.0 (- (* 4.0 (+ (count group1) (count group2))) 9.0)))]
    (* j (hedges-g group1 group2))))

(defn ameasure
  "Vargha-Delaney A measure for two populations a and b"
  {:metadoc/categories #{:effect}}
  ^double [group1 group2]
  (let [m (count group1)
        n (count group2)
        ^double r1 (reduce #(+ ^double %1 ^double %2) (take m (m/rank (concat group1 group2))))]
    (/ (- (+ r1 r1) (* m (inc m)))
       (* 2.0 m n))))

(defn cliffs-delta
  "Cliff's delta effect size"
  {:metadoc/categories #{:effect}}
  ^double [group1 group2]
  (/ ^double (reduce #(+ ^double %1 ^double %2)
                     (for [^double a group1
                           ^double b group2]
                       (m/signum (- a b))))
     (* (count group1) (count group2))))

;; binary classification statistics

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
    (let [f (cond
              (map? true-value) true-value
              (seqable? true-value) (partial contains? (set true-value))
              :else #(= % true-value))]
      (map f xs))))

(defn binary-measures-all
  "Collection of binary measures.

  * `truth` - list of ground truth values
  * `prediction` - list of predicted values
  * `true-value` - optional, what is true in `truth` and `prediction`

  `true-value` can be one of:

  * `nil` - values are treating as booleans
  * any sequence - values from sequence will be treated as `true`
  * map - conversion will be done according to provided map (if there is no correspondin key, value is treated as `false`)

  https://en.wikipedia.org/wiki/Precision_and_recall"
  {:metadoc/categories #{:stat}}
  ([truth prediction] (binary-measures-all truth prediction nil))
  ([truth prediction true-value]
   (let [truth (binary-process-list truth true-value)
         prediction (binary-process-list prediction true-value)
         {:keys [^double tp ^double fp ^double fn ^double tn] :as details} (merge {:tp 0.0 :fn 0.0 :fp 0.0 :tn 0.0}
                                                                                 (frequencies (map binary-confusion truth prediction)))
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
         f-beta (clojure.core/fn [^double beta] (let [b2 (* beta beta)]
                                                  (* (inc b2) (/ (* ppv tpr)
                                                                 (+ ppv tpr)))))
         f1-score (f-beta 1.0)]
     (merge details {:cp cp
                     :cn cn
                     :pcp pcp
                     :pcn pcn
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
                     :ppv ppv
                     :precision ppv
                     :fdr (- 1.0 ppv)
                     :npv npv
                     :for (- 1.0 npv)
                     :lr+ lr+
                     :lr- lr-
                     :dor (/ lr+ lr-)
                     :f-measure f1-score
                     :f1-score f1-score
                     :f-beta f-beta
                     :mcc (/ (- (* tp tn) (* fp fn))
                             (m/sqrt (* (+ tp fp)
                                        (+ tp fn)
                                        (+ tn fp)
                                        (+ tn fn))))
                     :bm (dec (+ tpr tnr))
                     :mk (dec (+ ppv npv))}))))

(defn binary-measures
  "Subset of binary measures. See [[binary-measures-all]].

  Following keys are returned: `[:tp :tn :fp :fn :accuracy :fdr :f-measure :fall-out :precision :recall :sensitivity :specificity :prevalance]`"
  {:metadoc/categories #{:stat}}
  ([truth prediction] (binary-measures truth prediction nil))
  ([truth prediction true-value]
   (select-keys (binary-measures-all truth prediction true-value)
                [:tp :tn :fp :fn :accuracy :fdr :f-measure :fall-out :precision :recall :sensitivity :specificity :prevalance])))

;; tests

;; t-test, reimplementation of R version

(defn- ttest-two-sided
  [^double tstat ^double alpha ^double df]
  (let [d (r/distribution :t {:degrees-of-freedom df})
        p (* 2.0 ^double (r/cdf d (- (m/abs tstat))))
        ^double cint (r/icdf d (- 1.0 (* 0.5 alpha)))]
    {:p-value p
     :confidence-intervals [(- tstat cint) (+ tstat cint)]}))

(defn- ttest-less
  [^double tstat ^double alpha ^double df]
  (let [d (r/distribution :t {:degrees-of-freedom df})]
    {:p-value (r/cdf d tstat)
     :confidence-intervals [##-Inf (+ tstat ^double (r/icdf d (- 1.0 alpha)))]}))

(defn- ttest-greater
  [^double tstat ^double alpha ^double df]
  (let [d (r/distribution :t {:degrees-of-freedom df})]
    {:p-value (- 1.0 ^double (r/cdf d tstat))
     :confidence-intervals [(- tstat ^double (r/icdf d (- 1.0 alpha))) ##Inf]}))

(defn- ttest-sides-fn
  [sides]
  (case sides
    :one-sided-less ttest-less
    :one-sided ttest-less
    :one-sided-greater ttest-greater
    ttest-two-sided))

(defn- ttest-update-ci
  [^double mu ^double stderr [^double l ^double r]]
  [(+ mu (* l stderr))
   (+ mu (* r stderr))])

(defn ttest-one-sample
  "One-sample Student's t-test

  * `alpha` - significance level (default: `0.05`)
  * `sides` - one of: `:two-sided`, `:one-sided-less` (short: `:one-sided`) or `:one-sided-greater`
  * `mu` - mean (default: `0.0`)"
  {:metadoc/categories #{:test}}
  ([xs] (ttest-one-sample xs {}))
  ([xs {:keys [^double alpha sides ^double mu]
        :or {alpha 0.05 sides :two-sided mu 0.0}}]
   (let [axs (m/seq->double-array xs)
         n (alength axs)
         m (mean axs)
         v (variance axs)
         stderr (m/sqrt (/ v n))]
     (assert (> stderr (* 10.0 m/MACHINE-EPSILON (m/abs m))) "Constant data, can't perform test.")
     (let [df (dec n)
           tstat (/ (- m mu) stderr)
           pvals (-> ((ttest-sides-fn sides) tstat alpha df)
                     (update :confidence-intervals (partial ttest-update-ci mu stderr)))]
       (merge pvals {:estimated-mu m
                     :df df
                     :t tstat
                     :test-type sides})))))

(defn- ttest-equal-variances
  [^double nx ^double ny ^double vx ^double vy]
  (let [df (- (+ nx ny) 2.0)
        v (/ (+ (* vx (dec nx))
                (* vy (dec ny))) df)]
    [df (m/sqrt (* v (+ (/ 1.0 nx)
                        (/ 1.0 ny))))]))

(defn- ttest-not-equal-variances
  [^double nx ^double ny ^double vx ^double vy]
  (let [stderrx (m/sqrt (/ vx nx))
        stderry (m/sqrt (/ vy ny))
        stderr (m/hypot-sqrt stderrx stderry)
        df (/ (m/sq (m/sq stderr))
              (+ (/ (m/sq (m/sq stderrx)) (dec nx))
                 (/ (m/sq (m/sq stderry)) (dec ny))))]
    [df stderr]))

(defn ttest-two-samples
  "Two-sample Student's t-test

  * `alpha` - significance level (default: `0.05`)
  * `sides` - one of: `:two-sided`, `:one-sided-less` (short: `:one-sided`) or `:one-sided-greater`
  * `mu` - mean (default: `0.0`)
  * `paired?` - unpaired or paired test, boolean (default: `false`)
  * `equal-variances?` - unequal or equal variances, boolean (default: `false`)"
  {:metadoc/categories #{:test}}
  ([xs ys] (ttest-two-samples xs ys {}))
  ([xs ys {:keys [^double alpha sides ^double mu paired? equal-variances?]
           :or {alpha 0.05 sides :two-sided mu 0.0 paired? false equal-variances? false}
           :as params}]
   (let [nx (count xs)
         ny (count ys)]
     (assert (or (and equal-variances? (< 2 (+ nx ny)) (pos? nx) (pos? ny))
                 (and (not equal-variances?)
                      (> nx 1) (> ny 1))) "Not enough observations.")
     (when paired? (assert (== nx ny) "Lengths of xs and ys should be equal."))
     (if paired? (-> (ttest-one-sample (map (fn [^double x ^double y] (- x y)) xs ys) params)
                     (assoc :paired? true))
         (let [axs (m/seq->double-array xs)
               ays (m/seq->double-array ys)
               mx (mean axs)
               my (mean ays)
               vx (variance axs)
               vy (variance ays)
               [df ^double stderr] (if equal-variances?
                                     (ttest-equal-variances nx ny vx vy)
                                     (ttest-not-equal-variances nx ny vx vy))
               tstat (/ (- mx my mu) stderr)
               pvals (-> ((ttest-sides-fn sides) tstat alpha df)
                         (update :confidence-intervals (partial ttest-update-ci mu stderr)))]
           (merge pvals {:estimated-mu [mx my]
                         :df df
                         :t tstat
                         :test-type sides
                         :paired? false}))))))

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
  ([data lags] (acf-ci data lags 0.05))
  ([data ^long lags ^double alpha]
   (let [acf-data (acf data lags)
         ci (p-acf-ci-value data alpha)]
     {:ci ci
      :acf acf-data
      :cis (map (fn [^double r]
                  (* ci (m/sqrt (dec (+ r r))))) (reductions (fn [^double acc ^double s]
                                                               (+ acc (* s s))) acf-data))})))


(defn- estimate-acceleration
  "Estimates acceleration for BCA bootstrap confidence interval computation"
  ([avs m]  (let [influence (map #(- ^double % m) avs)
                  num       (sum (map m/cb influence))
                  denom     (* 6 (m/sqrt (m/cb (sum (map m/sq influence)))))]
              (/  num denom))))

(defn- accelerated-quantile [^double z0 ^double z1 ^double a]
  (let [num (+ z0 z1)
        denom (- 1 (* a num))]
    (+ z0 (/ (double num) denom))))
(defn- empirical-cdf [vs value]
  (/ (double  (count (filter  #(<= ^double % value) vs))) (count vs)))


(defn percentile-bca-extent
  "Return bias corrected percentile range and mean for bootstrap samples. Also accounts for variance
   variations throught the accelaration parameter.
  See https://projecteuclid.org/euclid.ss/1032280214

  `p` - calculates extent of bias corrected `p` and `100-p` (default: `p=2.5`)"
  {:metadoc/categories #{:extent}}
  ([vs] (percentile-bca-extent vs 2.5))
  ([vs ^double p] (percentile-bca-extent vs p (- 100 p)))
  ([vs p1 p2] (percentile-bca-extent vs p1 p2 :legacy))
  ([vs ^double p1 ^double p2 estimation-strategy]
   (let [avs (m/seq->double-array vs)
         m (mean avs)
         accel (estimate-acceleration avs m)
        ;;  icdf of the number of bootstrap samples <= the mean
         ^double z0 (r/icdf r/default-normal (empirical-cdf vs m))
         ^double z1 (r/icdf r/default-normal (/ p1 100))
         ^double z2 (r/icdf r/default-normal (/ p2 100))
         q1 (r/cdf r/default-normal (accelerated-quantile z0 z1 accel))
         q2 (r/cdf r/default-normal (accelerated-quantile z0 z2 accel))]

     [(percentile avs (* 100  q1) estimation-strategy)
      (percentile avs (* 100  q2) estimation-strategy)
      m])))

(defn percentile-bc-extent
  "Return bias corrected percentile range and mean for bootstrap samples.
  See https://projecteuclid.org/euclid.ss/1032280214

  `p` - calculates extent of bias corrected `p` and `100-p` (default: `p=2.5`)"
  {:metadoc/categories #{:extent}}
  ([vs] (percentile-bc-extent vs 2.5))
  ([vs ^double p] (percentile-bc-extent vs p (- 100.0 p)))
  ([vs p1 p2] (percentile-bc-extent vs p1 p2 :legacy))
  ([vs ^double p1 ^double p2 estimation-strategy]
   (let [avs (m/seq->double-array vs)
         m (mean avs)
        ;;  icdf of the number of bootstrap samples <= the mean
         ^double z0 (r/icdf r/default-normal (empirical-cdf vs m))
         ^double z1 (r/icdf r/default-normal (/ p1 100))
         ^double z2 (r/icdf r/default-normal (/ p2 100))
         q1 (r/cdf r/default-normal (accelerated-quantile z0 z1 0.0))
         q2 (r/cdf r/default-normal (accelerated-quantile z0 z2 0.0))]

     [(percentile avs (* 100  q1) estimation-strategy)
      (percentile avs (* 100  q2) estimation-strategy)
      m])))
