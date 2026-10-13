(ns fastmath.signal.denoise
  "Wavelet shrinkage (denoising).

  [[threshold]] estimates a threshold from wavelet coefficients, [[denoise]] applies a shrinkage rule to the detail coefficients: either to the bands returned by `fastmath.transform/dwt`, to flat coefficients, or to a signal transformed with a given wavelet.

  Shrinkage rules for a threshold `t` (the first three as `pywt.threshold`):

  * `:hard` - `x` when `|x| >= t`, otherwise 0
  * `:soft` - `sign(x) * max(|x| - t, 0)`
  * `:garrote` - `x - t^2/x` when `|x| > t`, otherwise 0
  * `:hyperbole` - `sign(x) * sqrt(x^2 - t^2)` when `|x| > t`, otherwise 0"
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.vector :as v]
            [fastmath.protocols.wavelets :as prot]
            [fastmath.transform :as t])
  (:import [fastmath.java Array]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defn- non-empty!
  [xs]
  (when (empty? xs) (throw (ex-info "Coefficients are empty." {})))
  xs)

;; https://www.diva-portal.org/smash/get/diva2:1003644/FULLTEXT01.pdf

(defn sure
  "Stein's Unbiased Risk Estimate threshold for coefficients with unit noise level, as MATLAB `thselect(x, 'rigrsure')`.

  With `s` the sorted squared absolute values, the risk of the `k`-th candidate is `(n - 2k + sum(s_1..s_k) + (n - k) s_k) / n` and the threshold is `sqrt(s_k)` of the smallest risk.

  Throws `ex-info` for empty input."
  ^double [xs]
  (let [s (double-array (sort (map m/sq (non-empty! xs))))
        n (alength s)]
    (loop [k 1 cumsum 0.0 best-risk ##Inf best 0.0]
      (if (m/> k n)
        (m/sqrt best)
        (let [sk (Array/aget s (m/long-dec k))
              cumsum (m/+ cumsum sk)
              risk (m// (m/+ (m/- n (m/* 2.0 k)) cumsum (m/* (m/- n k) sk)) n)]
          (if (m/< risk best-risk)
            (recur (m/long-inc k) cumsum risk sk)
            (recur (m/long-inc k) cumsum best-risk best)))))))

;; https://computing.llnl.gov/sites/default/files/jei2001.pdf

(defn minfdr-threshold
  "Threshold from the `top-n` quantile of two-sided normal p-values of the coefficients: the value whose p-value equals that quantile.

  The noise level is [[robust-stddev]] when `robust?` is `true`, the standard deviation otherwise; a zero noise level gives threshold 0.

  Throws `ex-info` for empty input."
  ^double [xs ^double top-n robust?]
  (let [sigma (if robust? (stats/robust-stddev xs :mad) (stats/stddev (non-empty! xs)))]
    (if (m/zero? sigma)
      0.0
      (let [pm (-> (map (fn [^double x]
                          (m/* 2.0 (r/ccdf r/default-normal (m// (m/abs x) sigma)))) xs)
                   (stats/quantile top-n))]
        (m/* sigma (double (r/icdf r/default-normal (m/- 1.0 (m/* 0.5 pm)))))))))

(defn top-N-threshold
  "Threshold keeping the `top-n` fraction of the largest absolute coefficients: the `1 - top-n` quantile of absolute values."
  ^double [xs ^double top-n]
  (stats/quantile (v/abs (non-empty! xs)) (m/- 1.0 top-n)))

(def ^:private thresholds
  [:visu :universal :sure :rigrsure :hybrid :heursure :avg :peaksavg :peaks-avg :topn :top-n :minfdr :minfdr-robust])

(defn- scaled
  "`f` applied to coefficients divided by their [[robust-stddev]], scaled back; 0 when the deviation is 0."
  ^double [coeffs f]
  (let [sigma (stats/robust-stddev coeffs :mad)]
    (if (m/zero? sigma)
      0.0
      (m/* sigma (double (f (v/div coeffs sigma)))))))

(defn threshold
  "Estimates a denoising threshold from wavelet coefficients `coeffs`.

  `thr` is a number (returned as is) or one of:
  
  * `:visu` (default) - `sigma * sqrt(2 ln n)`, `sigma` estimated by [[robust-stddev]], `n` the number of coefficients
  * `:universal` - as `:visu`, with the standard deviation
  * `:sure` or `:rigrsure` - [[sure]] of the coefficients divided by `sigma` ([[robust-stddev]]), multiplied by `sigma`
  * `:hybrid` or `:heursure` - MATLAB `heursure` on the coefficients divided by `sigma`, multiplied by `sigma`: `sqrt(2 ln n)` when the signal energy is small, otherwise the smaller of SURE and `sqrt(2 ln n)`
  * `:avg` - mean of absolute values
  * `:peaksavg` or `:peaks-avg` - mid point between the minimum and maximum absolute value
  * `:topn` or `:top-n` - keep the `top-n-ratio` fraction of largest coefficients ([[top-N-threshold]])
  * `:minfdr`, `:minfdr-robust` - [[minfdr-threshold]] with the standard deviation or [[robust-stddev]]
  
  `top-n-ratio` (default 0.25) is used by `:topn`, `:minfdr` and `:minfdr-robust`. All given coefficients are used; pass only detail coefficients to estimate from them.

  Throws `ex-info` for empty coefficients or an unknown `thr`."
  (^double [coeffs] (threshold coeffs :visu))
  (^double [coeffs thr] (threshold coeffs thr 0.25))
  (^double [coeffs thr ^double top-n-ratio]
   (if (number? thr)
     (double thr)
     (let [n (count (non-empty! coeffs))
           sqrt-2-log-n (m/sqrt (m/* 2.0 (m/log n)))]
       (case thr
         :visu (m/* (stats/robust-stddev coeffs :mad) sqrt-2-log-n)
         :universal (m/* (stats/stddev coeffs) sqrt-2-log-n)
         (:sure :rigrsure) (scaled coeffs sure)
         (:hybrid :heursure) (scaled coeffs
                                     (fn [xs]
                                       (let [eta (m// (m/- (v/dot xs xs) n) n)
                                             crit (m// (m/pow (m/log2 n) 1.5) (m/sqrt n))]
                                         (if (m/< eta crit)
                                           sqrt-2-log-n
                                           (m/min (sure xs) sqrt-2-log-n)))))
         :avg (stats/mean (v/abs coeffs))
         (:peaksavg :peaks-avg) (stats/mean (stats/extent (v/abs coeffs) false))
         (:topn :top-n) (top-N-threshold coeffs top-n-ratio)
         :minfdr (minfdr-threshold coeffs top-n-ratio false)
         :minfdr-robust (minfdr-threshold coeffs top-n-ratio true)
         (throw (ex-info "Unknown threshold." {:thr thr :allowed thresholds})))))))

(def ^:private shrinkage-methods [:hard :soft :garrote :hyperbole])

(defn- shrink-fn
  [method ^double lambda]
  (let [l2 (m/* lambda lambda)]
    (case method
      :hard (fn ^double [^double x] (if (m/< (m/abs x) lambda) 0.0 x))
      :soft (fn ^double [^double x] (m/* (m/signum x) (m/max (m/- (m/abs x) lambda) 0.0)))
      :garrote (fn ^double [^double x] (if (m/> (m/abs x) lambda) (m/- x (m// l2 x)) 0.0))
      :hyperbole (fn ^double [^double x] (if (m/> (m/abs x) lambda) (m/* (m/signum x) (m/sqrt (m/- (m/* x x) l2))) 0.0))
      (throw (ex-info "Unknown shrinkage method." {:method method :allowed shrinkage-methods})))))

(defn- bands?
  [coeffs]
  (sequential? (first coeffs)))

(defn- flat-lengths
  "Band lengths of flat `dwt` coefficients, from their metadata."
  [coeffs]
  (:lengths (:fastmath.transform/dwt (meta coeffs))))

(defn- split-bands
  [coeffs lengths]
  (when-not (== (count coeffs) (long (reduce + lengths)))
    (throw (ex-info "Flat coefficients do not match band lengths." {:length (count coeffs) :lengths lengths})))
  (loop [cs (vec coeffs) ls lengths acc []]
    (if-let [l (first ls)]
      (recur (subvec cs l) (rest ls) (conj acc (subvec cs 0 l)))
      acc)))

(defn denoise
  "Applies wavelet shrinkage to the detail coefficients.

  Arities:

  * `[coeffs]`, `[coeffs options]` - `coeffs` are bands from `fastmath.transform/dwt` (`[cA cD_n ... cD_1]`) or flat coefficients
  * `[wavelet xs options]` - denoises signal `xs`: `fastmath.transform/dwt` with `wavelet` (a name or wavelet object) and `:level` / `:mode` from options, shrinkage of the bands, `fastmath.transform/idwt`; returns the signal of the original length
  * `[transformer xs options]` - as above with a transform object from `fastmath.transform/transformer` and flat coefficients without band lengths (the flat rule below, `:skip` protects the approximation); options are passed to the transform, the reverse transform also gets `:n`

  Options:

  * `:method` - shrinkage rule `:hard` (default), `:soft`, `:garrote` or `:hyperbole`, see the namespace docstring
  * `:thr` - a number or a method of [[threshold]] (default `:sure`)
  * `:top-n-ratio` - for [[threshold]] (default 0.25)
  * `:high?` - estimate the threshold from the finest details only (default `true` for `:visu` and `:universal`, `false` otherwise); for bands the finest detail band `cD_1`, for flat coefficients the second half
  * `:skip` - flat coefficients only: number of leading coefficients left unchanged (default 0)

  Bands: the approximation band is never changed; one threshold, estimated from the finest band or from all detail bands, is applied to every detail band. Flat coefficients from `fastmath.transform/dwt` with `:flat? true` (metadata with band `:lengths`) are treated as bands. Other flat coefficients: the threshold is estimated from the second half (`:high?`) or from all coefficients and applied to all coefficients after `:skip`.

  Returns bands, or a vector of flat coefficients, with the metadata of `coeffs` (so `fastmath.transform/idwt` still works); the signal arities return a vector of doubles. Throws `ex-info` for empty coefficients or an unknown method or threshold."
  ([coeffs] (denoise coeffs nil))
  ([coeffs {:keys [method thr skip high? top-n-ratio]
            :or {method :hard thr :sure skip 0 top-n-ratio 0.25}}]
   (non-empty! coeffs)
   (let [high? (if (nil? high?) (boolean (#{:visu :universal} thr)) high?)
         lengths (when-not (bands? coeffs) (flat-lengths coeffs))]
     (cond
       lengths
       (with-meta (vec (mapcat identity (denoise (split-bands coeffs lengths) {:method method :thr thr :high? high? :top-n-ratio top-n-ratio})))
         (meta coeffs))

       (bands? coeffs)
       (let [details (rest coeffs)
             estimate (if (or high? (empty? details))
                        (last coeffs)
                        (mapcat identity details))
             f (shrink-fn method (threshold estimate thr top-n-ratio))]
         (with-meta (into [(vec (first coeffs))] (map (fn [band] (mapv f band)) details))
           (meta coeffs)))

       :else
       (let [n (count coeffs)
             cs (vec coeffs)
             skip (long skip)
             estimate (if high? (subvec cs (m/long-quot n 2)) cs)
             f (shrink-fn method (threshold estimate thr top-n-ratio))]
         (with-meta (into (subvec cs 0 (m/long-min skip n)) (map f (subvec cs (m/long-min skip n))))
           (meta coeffs))))))
  ([wavelet-or-transformer xs options]
   (let [n (count (non-empty! xs))]
     (if (satisfies? prot/TransformProto wavelet-or-transformer)
       (let [opts (assoc options :n n)]
         (vec (prot/reverse-1d wavelet-or-transformer
                               (double-array (denoise (vec (prot/forward-1d wavelet-or-transformer xs options)) options))
                               opts)))
       (t/idwt (denoise (t/dwt xs wavelet-or-transformer options) options) wavelet-or-transformer)))))
