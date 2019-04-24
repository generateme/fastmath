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
                        :corr "Correlation"
                        :extent "Extents"}}
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.random :as rnd])
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

(defn sem
  "Standard error of mean"
  {:metadoc/categories #{:stat}}
  ^double [vs]
  (let [s (stddev vs)]
    (/ s (m/sqrt (count vs)))))

(defmacro ^:private build-extent
  [nm mid ext]
  `(defn ~nm
     ~(str mid " +/- " ext)
     {:metadoc/categories #{:extent}}
     [vs#]
     (let [vs# (m/seq->double-array vs#)
           m# (~mid vs#)
           s# (~ext vs#)]
       [(- m# s#) (+ m# s#) m#])))

(build-extent stddev-extent mean stddev)
(build-extent mad-extent median median-absolute-deviation)
(build-extent sem-extent mean sem)

(defn percentile-extent
  "Return percentile range."
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

(defn ci
  "T-student based confidence interval for given data. Alpha value defaults to 0.98."
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
  Last parameter is statistical function used to measure, default to mean."
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
         buff (long-array bins)
         samples (count vs)
         samplesd (double samples)]
     
     (doseq [^double v vs] 
       (let [b (java.util.Arrays/binarySearch ^doubles search-array v)
             ^int pos (if (neg? b) (m/abs (+ b 2)) b)]
         (fastmath.java.Array/inc ^longs buff pos)))

     {:size bins
      :step step
      :samples samples
      :min mn
      :max mx
      :bins (map vector search-array buff)})))

;;

(defonce ^:private ^:const ^double gaussian-factor (/ (m/sqrt m/TWO_PI)))

(defn- kde-uniform  ^double [^double x] (if (<= (m/abs x) 1.0) 0.5 0.0))
(defn- kde-gaussian ^double [^double x] (* gaussian-factor (m/exp (* -0.5 x x))))
(defn- kde-triangular ^double [^double x] (let [absx (m/abs x)]
                                            (if (<= absx 1.0) (- 1.0 absx) 0.0)))
(defn- kde-epanechnikov ^double [^double x] (if (<= (m/abs x) 1.0) (* 0.75 (- 1.0 (* x x))) 0.0))
(defn- kde-quartic ^double [^double x] (if (<= (m/abs x) 1.0) (* 0.9375 (m/sq (- 1.0 (* x x)))) 0.0))
(defn- kde-triweight
  ^double [^double x]
  (if (<= (m/abs x) 1.0)
    (let [v (- 1.0 (* x x))]
      (* 1.09375 v v v)) 0.0))

(defn- kde-tricube
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (<= absx 1.0)
      (let [v (- 1.0 (* absx absx absx))]
        (* 0.8875 v v v)) 0.0)))

(defn- kde-cosine ^double [^double x] (if (<= (m/abs x) 1.0) (* m/QUARTER_PI (m/cos (* m/HALF_PI x))) 0.0))
(defn- kde-logistic ^double [^double x] (/ (+ 2.0 (m/exp x) (m/exp (- x)))))
(defn- kde-sigmoid ^double [^double x] (/ (/ 2.0 (+ (m/exp x) (m/exp (- x)))) m/PI))
(defn- kde-silverman
  ^double [^double x]
  (let [xx (/ (m/abs x) m/SQRT2)]
    (* 0.5 (m/exp (- xx)) (m/sin (+ xx m/QUARTER_PI)))))

;; experimental
(defn- kde-laplace ^double [^double x] (* 0.5 (m/exp (- (m/abs x)))))
(defn- kde-wigner ^double [^double x] (if (<= (m/abs x) 1.0) (/ (* 2.0 (m/sqrt (- 1.0 (* x x)))) m/PI) 0.0))
(defn- kde-cauchy ^double [^double x] (/ (* m/HALF_PI (inc (m/sq (/ x 0.5))))))

(defn- nrd
  ^double [data]
  (let [sd (stats/stddev data)
        iqr (stats/iqr data)]
    (* 1.06 (min sd (/ iqr 1.34)) (m/pow (alength ^doubles data) -0.2))))

(defn- kde
  "Return kernel density estimation function"
  ([data k] (kde data k nil))
  ([data k h]
   (let [data (let [a (m/seq->double-array data)] (java.util.Arrays/sort a) a)
         h (double (or h (nrd data)))
         hrev (/ h)
         span (* 6.0 h)
         factor (/ (* (alength data) h))]
     [(fn [^double x]
        (let [start (java.util.Arrays/binarySearch data (- x span))
              ^int start (if (neg? start) (dec (- start)) start)
              end (java.util.Arrays/binarySearch data (+ x span))
              ^int end (if (neg? end) (dec (- end)) end)
              ^doubles xs (java.util.Arrays/copyOfRange data start end)]
          (* factor (double (areduce xs i sum (double 0.0)
                                     (+ sum ^double (k (* hrev (- x ^double (aget xs i))))))))))
      factor h])))

(defonce ^:private kde-integral
  {:uniform 0.5
   :triangular m/TWO_THIRD
   :epanechnikov 0.6
   :quartic (/ 5.0 7.0)
   :triweight (/ 350.0 429.0)
   :tricube (/ 175.0 247.0)
   :gaussian (* 0.5 (/ m/SQRTPI))
   :cosine (* 0.0625 m/PI m/PI)
   :logistic m/SIXTH
   :sigmoid (/ 2.0 (* m/PI m/PI))
   :silverman (* 0.0625 3.0 m/SQRT2)})

(defmulti kernel-density (fn [k & _] k))

(defmacro ^:private make-kernel-density-fns
  [lst]
  `(do ~@(for [v (eval lst)
               :let [n (symbol (str "kde-" (name v)))]]
           `(defmethod kernel-density ~v
              ([k# vs#] (first (kde vs# ~n)))
              ([k# vs# h#] (first (kde vs# ~n h#)))
              ([k# vs# h# all?#] (let [kded# (kde vs# ~n h#)]
                                   (if all?# kded# (first kded#))))))))

(make-kernel-density-fns (concat (keys kde-integral) [:wigner :laplace :cauchy]))

(defmethod kernel-density :smile
  ([_ vs h] (if h
              (let [^KernelDensity k (KernelDensity. (m/seq->double-array vs) h)]
                (fn [x] (.p k x)))
              (kernel-density :smile vs)))
  ([_ vs] (let [^KernelDensity k (KernelDensity. (m/seq->double-array vs))]
            (fn [x] (.p k x)))))

(defmethod kernel-density :default [_ & r] (apply kernel-density :smile r))

(defonce ^:private gaussian01 (r/distribution :normal))

(defn kernel-density-ci
  "Create function which returns confidence intervals for given kde method.

  Check 6.1.5 http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xlghtmlnode33.html"
  ([method data] (kernel-density-ci method data nil))
  ([method data bandwidth] (kernel-density-ci method data bandwidth 0.05))
  ([method data bandwidth ^double alpha]
   (if (contains? kde-integral method)
     (let [^double za (r/icdf gaussian01 (- 1.0 (* 0.5 (or alpha 0.05))))
           [kde-f ^double factor] (kernel-density method data bandwidth true)]
       (fn [^double x]
         (let [^double fx (kde-f x)
               band (* za (m/sqrt (* factor ^double (kde-integral method) fx)))]
           [fx (- fx band) (+ fx band)])))
     (let [kde-f (kernel-density method data bandwidth)]
       (fn [x] (let [fx (kde-f x)]
                [fx fx fx]))))))

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
  {:metadoc/categories #{:stat}}
  ^double [group1 group2]
  (cohens-d-with-correct group1 group2 -2.0))

(defn cohens-d
  "Cohen's d effect size for two groups"
  {:metadoc/categories #{:stat}}
  ^double [group1 group2]
  (cohens-d-with-correct group1 group2 0.0))

(defn glass-delta
  "Glass's delta effect size for two groups"
  {:metadoc/categories #{:stat}}
  ^double [group1 group2]
  (let [group2 (m/seq->double-array group2)]
    (/ (- (mean group1) (mean group2)) (stddev group2))))

(defn hedges-g
  "Hedges's g effect size for two groups"
  {:metadoc/categories #{:stat}}
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
  {:metadoc/categories #{:stat}}
  ^double [group1 group2]
  (let [j (- 1.0 (/ 3.0 (- (* 4.0 (+ (count group1) (count group2))) 9.0)))]
    (* j (hedges-g group1 group2))))

(defn- rank
  [vs]
  (let [m (into {} (map (fn [[k v]]
                          [k (/ ^double (reduce #(+ ^double %1 ^double %2) (map (comp #(inc ^double %) first) v))
                                (count v))]) (group-by second (map-indexed vector (sort vs)))))]
    (map m vs)))

(defn ameasure
  "Vargha-Delaney A measure for two populations a and b"
  ^double [group1 group2]
  (let [m (count group1)
        n (count group2)
        ^double r1 (reduce #(+ ^double %1 ^double %2) (take m (rank (concat group1 group2))))]
    (/ (- (+ r1 r1) (* m (inc m)))
       (* 2.0 m n))))

(defn cliffs-delta
  "Cliff's delta effect size"
  ^double [group1 group2]
  (/ ^double (reduce #(+ ^double %1 ^double %2)
                     (for [^double a group1
                           ^double b group2]
                       (m/signum (- a b))))
     (* (count group1) (count group2))))

#_(let [a1 (repeatedly 100 rand)
        a2 (repeatedly 1000 #(+ (rand) (rand)))]
    [(cohens-d a1 a2)
     (cohens-d-orig a1 a2)
     (glass-delta a1 a2)
     (hedges-g a1 a2)
     (hedges-g* a1 a2)
     (ameasure a1 a2)
     (cliffs-delta a1 a2)])


;;;

;; binary classification statistics

(defn- binary-confusion
  [t p]
  (cond
    (and t p) :tp
    (and t (not p)) :fn
    (and (not t) p) :fp
    :else :tn))

(defn- binary-confusion-val
  [tv fv t p]
  (cond
    (and (= t tv) (= p tv)) :tp
    (and (= t tv) (= p fv)) :fn
    (and (= t fv) (= p tv)) :fp
    :else :tn))

(defn binary-measures-all
  "https://en.wikipedia.org/wiki/Precision_and_recall"
  ([truth prediction] (binary-measures-all truth prediction nil nil))
  ([truth prediction true-value false-value]
   (let [confusion (if (and (nil? true-value) (nil? false-value))
                     binary-confusion
                     (partial binary-confusion-val true-value false-value))
         {:keys [^double tp ^double fp ^double fn ^double tn] :as details} (merge {:tp 0.0 :fn 0.0 :fp 0.0 :tn 0.0}
                                                                                 (frequencies (map confusion truth prediction)))
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
  ([truth prediction] (binary-measures truth prediction nil nil))
  ([truth prediction true-value false-value]
   (select-keys (binary-measures-all truth prediction true-value false-value)
                [:tp :tn :fp :fn :accuracy :fdr :f-measure :fall-out :precision :recall :sensitivity :specificity :prevalance])))


(comment defn t-test
         ""
         [sample1 sample2]
         (TestUtils/t (m/seq->double-array sample1) (m/seq->double-array sample2)))

(comment t-test [30.02 29.99 30.11 29.97 30.01 29.99]
         [29.89 29.93 29.72 29.98 30.02 29.98])

