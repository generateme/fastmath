(ns fastmath.stats
  "Statistics functions.

  * Descriptive statistics for sequence.
  * Correlation / covariance of two sequences.
  * Outliers

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
  * `:SecMoment` - second central moment, use: [[second-moment]]

  Note: [[percentile]] and [[quartile]] can have 10 different interpolation strategies. See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.html)

  ### Correlation / Covariance / Divergence

  * [[covariance]]
  * [[correlation]]
  * [[pearson-correlation]]
  * [[spearman-correlation]]
  * [[kendall-correlation]]
  * [[kullback-leibler-divergence]]
  * [[jensen-shannon-divergence]]

  ### Other

  Normalize samples to have mean=0 and standard deviation = 1 with [[standardize]].

  [[histogram]] to count samples in evenly spaced ranges."
  {:metadoc/categories {:stat "Descriptive statistics"
                        :corr "Correlation"}}
  (:require [fastmath.core :as m])
  (:import [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.stat.descriptive.rank Percentile Percentile$EstimationType]
           [org.apache.commons.math3.stat.descriptive.moment Kurtosis SecondMoment Skewness]
           [org.apache.commons.math3.stat.correlation KendallsCorrelation SpearmansCorrelation PearsonsCorrelation]
           [smile.stat.distribution KernelDensity]
           [org.apache.commons.math3.stat.inference TestUtils]))

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

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/index.html)
  
  See also [[quantile]]."
  {:metadoc/categories #{:stat}}
  (^double [vs ^double p] 
   (StatUtils/percentile (m/seq->double-array vs) p))
  (^double [vs ^double p estimation-strategy]
   (let [^Percentile perc (.withEstimationType ^Percentile (Percentile.) (or (estimation-strategies-list estimation-strategy) Percentile$EstimationType/LEGACY))]
     (.evaluate perc (m/seq->double-array vs) p ))))

(defn quantile
  "Calculate quantile of a `vs`.

  Percentile `p` is from range 0.0-1.0.
  
  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html) for interpolation strategy.

  Optionally you can provide `estimation-strategy` to change interpolation methods for selecting values. Default is `:legacy`. See more [here](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/index.html)
  
  See also [[percentile]]."
  {:metadoc/categories #{:stat}}
  (^double [vs ^double p]
   (percentile vs (m/constrain (* p 100.0) 0.0 100.0)))
  (^double [vs ^double p estimation-strategy]
   (percentile vs (m/constrain (* p 100.0) 0.0 100.0) estimation-strategy)))

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
  (smile.math.Math/mad (m/seq->double-array vs)))

(defn iqr
  "Interquartile range."
  {:metadoc/categories #{:stat}}
  (^double [vs] (iqr vs :legacy))
  (^double [vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         q1 (percentile avs 25.0 estimation-strategy)
         q3 (percentile avs 75.0 estimation-strategy)]
     (- q3 q1))))

(defn adjacent-values
  "Lower and upper adjacent values (LAV and UAV).

  Let Q1 is 25-percentile and Q3 is 75-percentile. IQR is `(- Q3 Q1)`.

  * LAV is smallest value which is greater or equal to the LIF = `(- Q1 (* 1.5 IQR))`.
  * UAV is largest value which is lower or equal to the UIF = `(+ Q3 (* 1.5 IQR))`.

  Optional `estimation-strategy` argument can be set to change quantile calculations estimation type. See [[estimation-strategies]]."
  {:metadoc/categories #{:stat}}
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
      (last (filter #(<= (double %) uav-thr) avs))])))

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
   (let [avs (m/seq->double-array vs)
         iqr (* 1.5 (- q3 q1))
         lof-thr (- q1 iqr)
         uof-thr (+ q3 iqr)]
     (java.util.Arrays/sort avs)
     (filter #(let [v (double %)]
                (bool-or (< v lof-thr)
                         (> v uof-thr))) avs))))

(defn minimum
  "Minimum value from sequence."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (smile.math.Math/min ^doubles vs)
    (reduce clojure.core/min vs)))

(defn maximum
  "Maximum value from sequence."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (smile.math.Math/max ^doubles vs)
    (reduce clojure.core/max vs)))

(defn extent
  "Return extent (min, max) values from sequence"
  {:metadoc/categories #{:stat}}
  [vs]
  (let [^double fv (first vs)]
    (reduce (fn [[^double mn ^double mx] ^double v]
              [(min mn v) (max mx v)]) [fv fv] (rest vs))))

(defn sum
  "Sum of all `vs` values."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (smile.math.Math/sum ^doubles vs)
    (reduce clojure.core/+ vs)))

(defn kurtosis
  "Calculate kurtosis from sequence."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [^Kurtosis k (Kurtosis.)]
    (.evaluate k (m/seq->double-array vs))))

(defn second-moment
  "Calculate second moment from sequence.

  It's a sum of squared deviations from the sample mean"
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [^SecondMoment k (SecondMoment.)]
    (.evaluate k (m/seq->double-array vs))))

(defn skewness
  "Calculate kurtosis from sequence."
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [^Skewness k (Skewness.)]
    (.evaluate k (m/seq->double-array vs))))

(defn stats-map
  "Calculate several statistics of `vs` and return as map.

  Optional `estimation-strategy` argument can be set to change quantile calculations estimation type. See [[estimation-strategies]]."
  {:metadoc/categories #{:stat}}
  ([vs] (stats-map vs :legacy))
  ([vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         sz (alength avs)
         mn (smile.math.Math/min avs)
         mx (smile.math.Math/max avs)
         sm (smile.math.Math/sum avs)
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
      :Skewness (skewness avs)
      :SecMoment (second-moment avs)})))

(defn standardize
  "Normalize samples to have mean = 0 and stddev = 1."
  [vs]
  (seq ^doubles (StatUtils/normalize (m/seq->double-array vs))))

(defn covariance
  "Covariance of two sequences."
  {:metadoc/categories #{:corr}}
  [vs1 vs2]
  (smile.math.Math/cov (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn correlation
  "Correlation of two sequences."
  {:metadoc/categories #{:corr}}
  [vs1 vs2]
  (smile.math.Math/cor (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn spearman-correlation
  "Spearman's correlation of two sequences."
  {:metadoc/categories #{:corr}}
  [vs1 vs2]
  (.correlation ^SpearmansCorrelation (SpearmansCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn pearson-correlation
  "Pearson's correlation of two sequences."
  {:metadoc/categories #{:corr}}
  [vs1 vs2]
  (.correlation ^PearsonsCorrelation (PearsonsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn kendall-correlation
  "Kendall's correlation of two sequences."
  {:metadoc/categories #{:corr}}
  [vs1 vs2]
  (.correlation ^KendallsCorrelation (KendallsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn kullback-leibler-divergence
  "Kullback-Leibler divergence of two sequences."
  {:metadoc/categories #{:corr}}
  [vs1 vs2]
  (smile.math.Math/KullbackLeiblerDivergence (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn jensen-shannon-divergence
  "Jensen-Shannon divergence of two sequences."
  {:metadoc/categories #{:corr}}
  [vs1 vs2]
  (smile.math.Math/JensenShannonDivergence (m/seq->double-array vs1) (m/seq->double-array vs2)))

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
  ([vs method]
   (if-not (keyword? method)
     (or method (estimate-bins vs))
     (let [n (count vs)]
       (int (condp = method
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
              :freedman-diaconis (let [vvs (m/seq->double-array vs)
                                       h (/ (* 2.0 (iqr vvs))
                                            (m/cbrt n))]
                                   (scott-fd-helper vvs h))))))))

(defn histogram
  "Calculate histogram.

  Returns map with keys:

  * `:size` - number of bins
  * `:step` - distance between bins
  * `:bins` - list of triples of range lower value, number of hits and ratio of used samples
  * `:min` - min value
  * `:max` - max value
  * `:samples` - number of used samples

  For estimation methods check [[estimate-bins]]."
  {:metadoc/categories #{:stat}}
  ([vs] (histogram vs :freedman-diaconis))
  ([vs bins-or-estimate-method] (histogram vs (estimate-bins vs bins-or-estimate-method) (extent vs)))
  ([vs ^long bins [^double mn ^double mx]]
   (let [diff (- mx mn)
         step (/ diff bins)
         search-array (double-array (map #(+ mn (* ^long % step)) (range bins)))
         buff (long-array bins)
         mx+ (m/next-double mx)
         vs- (filter #(<= mn ^double % mx+) vs)
         samples (count vs-)
         samplesd (double samples)]
     
     (doseq [^double v vs-] 
       (let [b (java.util.Arrays/binarySearch ^doubles search-array v)
             ^int pos (if (neg? b) (m/abs (+ b 2)) b)]
         (fastmath.java.Array/inc ^longs buff pos)))

     {:size bins
      :step step
      :samples samples
      :min mn
      :max mx
      :bins (map #(vector %1 %2 (/ ^long %2 samplesd)) search-array buff)})))

;;

(defn kernel-density
  "Creates kernel density function for given series `vs` and optional bandwidth `h`."
  {:metadoc/categories #{:stat}}
  ([vs ^double h]
   (let [^KernelDensity k (KernelDensity. (m/seq->double-array vs) h)]
     (fn [x] (.p k x))))
  ([vs]
   (let [^KernelDensity k (KernelDensity. (m/seq->double-array vs))]
     (fn [x] (.p k x)))))

;;;;;;;;;;;;;;
;; tests

(comment defn t-test
         ""
         [sample1 sample2]
         (TestUtils/t (m/seq->double-array sample1) (m/seq->double-array sample2)))

(comment t-test [30.02 29.99 30.11 29.97 30.01 29.99]
         [29.89 29.93 29.72 29.98 30.02 29.98])
