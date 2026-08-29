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
  "Calculate number of bins based on width of the bin."
  ^long [vvs ^double h]
  (let [h (if (m/<= 0.0 h m/EPSILON) (mad vvs) h)
        fv (first vvs)
        ^Vec2 mm (reduce (fn [^Vec2 curr ^double v]
                           (Vec2. (m/min (.x curr) v) (m/max (.y curr) v))) (Vec2. fv fv) (rest vvs))]
    (m/max 1 (long (m/ceil (m// (m/- (.y mm) (.x mm)) h))))))

(defn sturges
  ^long [^long n]
  (m/max 1 (long (m/inc (m/ceil (m/log2 n))))))

(defn rice
  ^long [^long n]
  (m/max 1 (long (m/ceil (m/* 2.0 (m/cbrt n))))))

(defn doane
  ^long [^doubles avs ^long n]
  (if (m/< n 3)
    1
    (let [stats (DescriptiveStatistics. avs)]
      (m/max 1 (m/+ (m/inc (m/log2 n))
                    (m/log2 (m/inc (m// (m/abs (.getSkewness stats))
                                        (m/sqrt (m// (m/* 6.0 (m/- n 2.0))
                                                     (m/* (m/inc n) (m/+ n 3.0))))))))))))

(defn scott
  ^long [^doubles avs ^long n]
  (let [h (m// (m/* 3.5 (m/sqrt (StatUtils/variance avs)))
               (m/cbrt n))]
    (scott-fd-helper avs h)))

(defn freedman-diaconis
  ^long [^doubles avs ^long n]
  (let [h (m// (m/* 2.0 (m/- (StatUtils/percentile avs 75)
                             (StatUtils/percentile avs 25)))
               (m/cbrt n))]
    (scott-fd-helper avs h)))
