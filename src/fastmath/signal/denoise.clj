(ns fastmath.signal.denoise
  "Wavelet denoise and shrinkage methods."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.vector :as v]
            [fastmath.protocols.wavelets :as prot])
  (:import [fastmath.java Array]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; https://www.diva-portal.org/smash/get/diva2:1003644/FULLTEXT01.pdf

(defn sure
  "Calculates SURE thresholding value."
  (^double [xs] (sure xs (count xs)))
  (^double [xs ^long n]
   (let [sxs (double-array (sort (map m/abs xs)))]
     (loop [k (long 1)
            curr (Array/get sxs 0)
            s (m/* curr curr)
            minrisk ##Inf]
       (if (m/> k n)
         curr
         (let [v (Array/get sxs (dec k))
               v2 (m/* v v)
               risk (m/+ (m/- n (m/* 2 k)) s (m/* (m/- n k) v2))]
           (if (m/< risk minrisk)
             (recur (inc k) v (m/+ s v2) risk)
             (recur (inc k) curr (m/+ s v2) minrisk))))))))

;; https://computing.llnl.gov/sites/default/files/jei2001.pdf

(defn robust-stddev
  "Stddev estimation based on MAD."
  ^double [xs]
  (-> (stats/median-absolute-deviation xs)
      (m// 0.6744897501960816))) ;; icdf(normal,0.75)

(defn minfdr-threshold
  "Minimal false discovery rate."
  ^double [xs ^double top-n robust?]
  (let [sigma (if robust? (robust-stddev xs) (stats/stddev xs))
        pm (-> (map (fn [^double x]
                      (m/* 2.0 (r/ccdf r/default-normal (m// (m/abs x) sigma)))) xs)
               (stats/quantile top-n))]
    (m/* sigma (double (r/icdf r/default-normal (m/- 1.0 (m/* 0.5 pm)))))))

(defn top-N-threshold
  "Keep Top N."
  ^double [xs ^double top-n]
  (stats/quantile (v/abs xs) (m/- 1.0 top-n)))

(defn threshold
  "Calculate optimal denoise threshold for wavelet coefficients.

  `thr` is one of the following
  
  * `:visu` - based on median absolute deviation estimate (default)
  * `:universal` - based on standard deviation estimate
  * `:sure` or `:rigrsure` - based on SURE estimator
  * `:hybrid` or `:heursure` - hybrid SURE estimator
  * `:avg` - abs coefficients average
  * `:peaksavg` - mid point between min and max absolute value of coefficents
  * `:topn` - keep top n/N largest coefficients
  * `:minfdr`, `minfdr-robust` - keep top n/N coefficients based on p-values (robust - uses MAD for sigma estimation)
  
  `top-n-ratio` determines cut-off point from `:topn`, `:minfdr` and `:midfdr-robust` thresholds.

  Always uses all coefficients data, remove first half of coefficients to use only details level."
  (^double [coeffs] (threshold coeffs :visu))
  (^double [coeffs thr] (threshold coeffs thr 0.25))
  (^double [coeffs thr ^double top-n-ratio]
   (let [n (count coeffs)]
     (if (number? thr)
       thr
       (case thr
         :visu (-> (robust-stddev coeffs)
                   (m/* (m/sqrt (m/* 2.0 (m/log n)))))
         :universal (-> (stats/stddev coeffs)
                        (m/* (m/sqrt (m/* 2.0 (m/log n)))))
         (:sure :rigrsure) (sure coeffs n)
         (:hybrid :heursure) (let [eta (m// (m/- (v/dot coeffs coeffs) n) n)
                                   crit (m// (m/pow (m/log2 n) 1.5) (m/sqrt n))]
                               (if (m/< eta crit)
                                 (m/sqrt (m/* 2.0 (m/log n)))
                                 (m/min (sure coeffs n) (m/sqrt (m/* 2.0 (m/log n))))))
         :avg (stats/mean (v/abs coeffs))
         (:peaksavg :peaks-avg) (stats/mean (stats/extent (v/abs coeffs) false))
         (:topn :top-n) (top-N-threshold coeffs top-n-ratio)
         :minfdr (minfdr-threshold coeffs top-n-ratio false)
         :minfdr-robust (minfdr-threshold coeffs top-n-ratio true))))))

(defn denoise
  "Wavelet shrinkage with some threshold for 1d signals.

  `coeffs` is a result of `dwt` or `wpt` transforms.

  Options:  

  * `:method` can be one of the following:
      * `:hard` (default)  
      * `:soft`
      * `:garrote`
      * `:hyperbole`
  * `:thr` can be a number of one of the [[threshold]] methods (default: `:sure`)
  * `:skip` can be used to leave `:skip` number of coefficients unaffected (default: 0)
  * `:high?` use only details level to estimate a threshold (`true` for `:visu` and `:universal`, `false` otherwise)
  * `:top-n-ratio` for `:topn` and `minfdr` sets ratio of highest coefficients to keep  

  Use on transformed sequences or call with transformer object."
  ([coeffs {:keys [method thr ^long skip high? ^double top-n-ratio]
            :or {method :hard thr :sure skip 0 top-n-ratio 0.25}}]
   (let [high? (if (nil? high?) (#{:visu :universal} thr) high?)
         n (count coeffs)
         t (double-array coeffs)
         lambda (threshold (if high? (drop (m// n 2) coeffs) coeffs) thr top-n-ratio)
         ids (range skip n)]
     (case method
       :soft (doseq [^long i ids]
               (let [v (Array/aget t i)]
                 (Array/aset t i (m/* (m/signum v) (m/max (m/- (m/abs v) lambda) 0.0)))))
       :hard (doseq [^long i ids]
               (let [v (Array/aget t i)]
                 (when (m/< (m/abs v) lambda) (Array/aset t i 0.0))))
       :garrote (let [l2 (m/sq lambda)]
                  (doseq [^long i ids]
                    (let [v (Array/aget t i)]
                      (Array/aset t i (if (m/> (m/abs v) lambda)
                                        (m/- v (m// l2 v))
                                        0.0)))))
       :hyperbole (let [l2 (m/sq lambda)]
                    (doseq [^long i ids]
                      (let [v (Array/aget t i)]
                        (Array/aset t i (if (m/> (m/abs v) lambda)
                                          (m/* (m/signum v) (m/sqrt (m/- (m/* v v) l2)))
                                          0.0))))))
     t))
  ([trans xs method]
   (let [v (prot/forward-1d trans xs)]
     (prot/reverse-1d trans (denoise v method))))
  ([coeffs] (denoise coeffs 0)))
