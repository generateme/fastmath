(ns fastmath.stats
  "Statistics functions.

  * Descriptive statistics for sequence.
  * Correlation / covariance of two sequences.
  * Outliers

  All functions are backed by Apache Commons Math or SMILE libraries. All work with Clojure sequences."
  {:metadoc/categories {:stat "Descriptive statistics"
                        :corr "Correlation"}}
  (:require [metadoc.examples :refer :all]
            [fastmath.core :as m]
            [fastmath.random :as r])
  (:import [org.apache.commons.math3.stat StatUtils]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn mode
  "Find the value that appears most often in a dataset `vs`.

  See also [[modes]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Example" (mode [1 2 3 -1 -1 2 -1 11 111]))
                      (example "Returns lowest value when every element appears equally." (mode [5 1 2 3 4]))]}
  ^double [vs]
  (let [m (StatUtils/mode (m/seq->double-array vs))]
    (aget ^doubles m 0)))

(defn modes
  "Find the values that appears most often in a dataset `vs`.

  Returns sequence with all most appearing values in increasing order.

  See also [[mode]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Example" (modes [1 2 3 -1 -1 2 -1 11 111]))
                      (example "Returns lowest value when every element appears equally." (modes [5 5 1 1 2 3 4 4]))]}
  [vs]
  (seq ^doubles (StatUtils/mode (m/seq->double-array vs))))

(defn percentile
  "Calculate percentile of a `vs`.

  Percentile `p` is from range 0-100.
  
  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html) for interpolation strategy.

  See also [[quantile]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Percentile 25%" (percentile [1 2 3 -1 -1 2 -1 11 111] 25.0))
                      (example "Percentile 50% (median)" (percentile [1 2 3 -1 -1 2 -1 11 111] 50.0))
                      (example "Percentile 75%" (percentile [1 2 3 -1 -1 2 -1 11 111] 75.0))
                      (example "Percentile 90%" (percentile [1 2 3 -1 -1 2 -1 11 111] 90.0))]}
  ^double [vs ^double p] 
  (StatUtils/percentile (m/seq->double-array vs) p))

(defn quantile
  "Calculate quantile of a `vs`.

  Percentile `p` is from range 0.0-1.0.
  
  See [docs](http://commons.apache.org/proper/commons-math/javadocs/api-3.4/org/apache/commons/math3/stat/descriptive/rank/Percentile.html) for interpolation strategy.

  See also [[percentile]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Quantile 0.25" (quantile [1 2 3 -1 -1 2 -1 11 111] 0.25))
                      (example "Quantile 0.5 (median)" (quantile [1 2 3 -1 -1 2 -1 11 111] 0.5))
                      (example "Quantile 0.75" (quantile [1 2 3 -1 -1 2 -1 11 111] 0.75))
                      (example "Quantile 0.9" (quantile [1 2 3 -1 -1 2 -1 11 111] 0.9))]}
  ^double [vs ^double p]
  (percentile vs (m/constrain (* p 100.0) 0.0 100.0)))

(defn median
  "Calculate median of a list. See [[median-3]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Median (percentile 50%)." (median [1 2 3 -1 -1 2 -1 11 111]))
                      (example "For three elements use faster [[median-3]]." (median [7 1 4]))]}
  ^double [vs]
  (percentile vs 50.0))

(defn median-3
  "Median of three values. See [[median]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Median of [7 1 4]" (median-3 7 1 4))]}
  ^double [^double a ^double b ^double c]
  (m/max (m/min a b) (m/min (m/max a b) c)))

(defn mean
  "Calculate mean of a list"
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Mean (average value)" (mean [1 2 3 -1 -1 2 -1 11 111]))]}
  (^double [vs] (StatUtils/mean (m/seq->double-array vs))))

(defn population-variance
  "Calculate population variance of a list.

  See [[variance]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Population variance" (population-variance [1 2 3 -1 -1 2 -1 11 111]))]}
  (^double [vs]
   (StatUtils/populationVariance (m/seq->double-array vs)))
  (^double [vs ^double u]
   (StatUtils/populationVariance (m/seq->double-array vs) u)))

(defn variance
  "Calculate variance of a list.

  See [[population-variance]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Variance." (variance [1 2 3 -1 -1 2 -1 11 111]))]}
  (^double [vs]
   (StatUtils/variance (m/seq->double-array vs)))
  (^double [vs ^double u]
   (StatUtils/variance (m/seq->double-array vs) u)))

(defn population-stddev
  "Calculate population standard deviation of a list.

  See [[stddev]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Population standard deviation." (population-stddev [1 2 3 -1 -1 2 -1 11 111]))]}
  (^double [vs]
   (m/sqrt (population-variance vs)))
  (^double [vs u]
   (m/sqrt (population-variance vs u))))

(defn stddev
  "Calculate population standard deviation of a list.

  See [[population-stddev]]."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Standard deviation." (stddev [1 2 3 -1 -1 2 -1 11 111]))]}
  (^double [vs]
   (m/sqrt (variance vs)))
  (^double [vs u]
   (m/sqrt (variance vs u))))

(defn median-absolute-deviation 
  "Calculate MAD"
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "MAD" (median-absolute-deviation [1 2 3 -1 -1 2 -1 11 111]))]}
  ^double [vs]
  (smile.math.Math/mad (m/seq->double-array vs)))

(defn adjacent-values
  "Lower and upper adjacent values (LAV and UAV).

  Let Q1 is 25-percentile and Q3 is 75-percentile. IQR is `(- Q3 Q1)`.

  * LAV is smallest value which is greater or equal to the LIF = `(- Q1 (* 1.5 IQR))`.
  * UAV is largest value which is lower or equal to the UIF = `(+ Q3 (* 1.5 IQR))`."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "[LAV, UAV]" (adjacent-values [1 2 3 -1 -1 2 -1 11 111]))
                      (example "Gaussian distribution [LAV, UAV]" (adjacent-values (repeatedly 1000000 r/grand)))]}
  ([vs]
   (let [avs (m/seq->double-array vs)
         q1 (percentile avs 25.0)
         q3 (percentile avs 75.0)]
     (adjacent-values avs q1 q3)))
  ([vs ^double q1 ^double q3]
   (let [avs (double-array vs)
         iqr (* 1.5 (- q3 q1))
         lav-thr (- q1 iqr)
         uav-thr (+ q3 iqr)]
     (java.util.Arrays/sort avs)
     [(first (filter #(>= (double %) lav-thr) avs))
      (last (filter #(<= (double %) uav-thr) avs))])))

(defn outliers
  "Find outliers defined as values outside outer fences.

  Let Q1 is 25-percentile and Q3 is 75-percentile. IQR is `(- Q3 Q1)`.

  * LOF (Lower Outer Fence) equals `(- Q1 (* 3.0 IQR))`.
  * UOF (Upper Outer Fence) equals `(+ Q3 (* 3.0 IQR))`.

  Returns sequence."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Outliers" (outliers [1 2 3 -1 -1 2 -1 11 111]))
                      (example "Gaussian distribution outliers" (outliers (repeatedly 3000000 r/grand)))]}
  ([vs]
   (let [avs (m/seq->double-array vs)
         q1 (percentile avs 25.0)
         q3 (percentile avs 75.0)]
     (outliers avs q1 q3)))
  ([vs ^double q1 ^double q3]
   (let [avs (double-array vs)
         iqr (* 3.0 (- q3 q1))
         lof-thr (- q1 iqr)
         uof-thr (+ q3 iqr)]
     (java.util.Arrays/sort avs)
     (filter #(let [v (double %)]
                (bool-or (< v lof-thr)
                         (> v uof-thr))) avs))))

(defn minimum
  "Minimum value from sequence."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Minimum value" (minimum [1 2 3 -1 -1 2 -1 11 111]))]}
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (smile.math.Math/min ^doubles vs)
    (reduce clojure.core/min vs)))

(defn maximum
  "Maximum value from sequence."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Maximum value" (maximum [1 2 3 -1 -1 2 -1 11 111]))]}
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (smile.math.Math/max ^doubles vs)
    (reduce clojure.core/max vs)))

(defn sum
  "Sum of all `vs` values."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Sum" (sum [1 2 3 -1 -1 2 -1 11 111]))]}
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (smile.math.Math/sum ^doubles vs)
    (reduce clojure.core/+ vs)))

(defn stats-map
  "Calculate several statistics from the list and return as map."
  {:metadoc/categories #{:stat}
   :metadoc/examples [(example "Stats" (stats-map [1 2 3 -1 -1 2 -1 11 111]))
                      (example "Select keys" (stats-map [:Mean :Q1 :Q3] [1 2 3 -1 -1 2 -1 11 111]))]}
  ([vs]
   (let [avs (m/seq->double-array vs)
         sz (alength avs)
         mn (smile.math.Math/min avs)
         mx (smile.math.Math/max avs)
         sm (smile.math.Math/sum avs)
         u (/ sm sz)
         mdn (median avs)
         q1 (percentile avs 25.0)
         q3 (percentile avs 75.0)
         sd (population-stddev avs)
         mad (median-absolute-deviation avs)
         [lav uav] (adjacent-values avs q1 q3)
         outliers (count (outliers avs q1 q3))]
     {:Size sz
      :Min mn
      :Max mx
      :Mean u
      :Median mdn
      :Mode (mode avs)
      :Q1 q1
      :Q3 q3
      :Total sm
      :SD sd
      :MAD mad
      :LAV lav
      :UAV uav
      :Outliers outliers}))
  ([ks vs]
   (zipmap ks (map (stats-map vs) ks))))

(defn standardize
  "Normalize samples to have mean = 0 and stddev = 1."
  {:metadoc/examples [(example "Standardize" (standardize [1 2 3 -1 -1 2 -1 11 111]))]}
  [vs]
  (seq ^doubles (StatUtils/normalize (m/seq->double-array vs))))

(defn covariance
  "Covariance of two sequences."
  {:metadoc/categories #{:corr}
   :metadoc/examples [(example "Covariance of uniform and gaussian distribution samples." (covariance (repeatedly 100000 r/grand) (repeatedly 100000 r/drand)))]}
  [vs1 vs2]
  (smile.math.Math/cov (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn correlation
  "Correlation of two sequences."
  {:metadoc/categories #{:corr}
   :metadoc/examples [(example "Correlation of uniform and gaussian distribution samples." (correlation (repeatedly 100000 r/grand) (repeatedly 100000 r/drand)))]}
  [vs1 vs2]
  (smile.math.Math/cor (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn spearman-correlation
  "Spearman correlation of two sequences."
  {:metadoc/categories #{:corr}
   :metadoc/examples [(example "Spearman correclation of uniform and gaussian distribution samples." (spearman-correlation (repeatedly 100000 r/grand) (repeatedly 100000 r/drand)))]}
  [vs1 vs2]
  (smile.math.Math/spearman (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn kendall-correlation
  "Kendall correlation of two sequences. Slow."
  {:metadoc/categories #{:corr}
   :metadoc/examples [(example "Kendall correlation of uniform and gaussian distribution samples." (kendall-correlation (repeatedly 100 r/grand) (repeatedly 100 r/drand)))]}
  [vs1 vs2]
  (smile.math.Math/kendall (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn kullback-leibler-divergence
  "Kullback-Leibler divergence of two sequences."
  {:metadoc/categories #{:corr}
   :metadoc/examples [(example "Kullback-Leibler divergence." (kullback-leibler-divergence (repeatedly 100 #(r/irand 100)) (repeatedly 100 #(r/irand 100))))]}
  [vs1 vs2]
  (smile.math.Math/KullbackLeiblerDivergence (m/seq->double-array vs1) (m/seq->double-array vs2)))

(defn jensen-shannon-divergence
  "Jensen-Shannon divergence of two sequences."
  {:metadoc/categories #{:corr}
   :metadoc/examples [(example "Jensen-Shannon divergence" (jensen-shannon-divergence (repeatedly 100 #(r/irand 100)) (repeatedly 100 #(r/irand 100))))]}
  [vs1 vs2]
  (smile.math.Math/JensenShannonDivergence (m/seq->double-array vs1) (m/seq->double-array vs2)))

;; TODO - replace with native (SMILE or Apache Commens) algorithms

(defn- closest-mean-fn
  [means]
  (fn [^double v] (reduce (partial min-key #(m/sq (- v ^double %))) means)))

;; `(k-means 4 '(1 2 3 -1 -1 2 -1 11 111)) => (-1.0 2.0 11.0 111.0)`
(defn k-means
  "k-means clustering"
  {:metadoc/examples [(example "Reduce to 4 values." (k-means 4 [1 2 3 -1 -1 2 -1 11 111]))]}
  [^long k vs]
  (let [vs (map double vs)
        svs (set vs)]
    (if (> k (count svs))
      (sort svs)
      (loop [mns (sort (take k (shuffle svs)))
             pmns (repeat k Double/NaN)]
        (if (= mns pmns)
          mns
          (recur (sort (map mean (vals (group-by (closest-mean-fn mns) vs)))) mns))))))
