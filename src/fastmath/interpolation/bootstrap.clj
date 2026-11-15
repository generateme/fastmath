(ns fastmath.interpolation.bootstrap
  "Bootstrap bands of one dimensional interpolators and smoothers: the fit refitted on bootstrap samples of the data
  pairs (drawn by [[fastmath.stats.bootstrap/bootstrap]]); the band at a point - a bootstrap confidence interval
  ([[fastmath.stats.bootstrap]] `ci-*`) of the refitted values there."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.interpolation.common :as ic])
  (:import [clojure.lang Fn IFn]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:private true :const true :tag 'long} default-samples 200)

(def intervals
  "Bootstrap confidence intervals of a band point: `:ci` keyword -> `(fn [boot-data alpha estimation] [lo hi t0])`
  over the fit at the point (`:t0`) and the refitted values (`:ts`), see [[fastmath.stats.bootstrap]]. Intervals
  needing more than these (`ci-bca`, `ci-studentized`, `ci-t`) don't apply."
  {:percentile boot/ci-percentile
   :basic boot/ci-basic
   :normal (fn [boot-data alpha _] (boot/ci-normal boot-data alpha))
   :bc boot/ci-bc})

(defn- check-level!
  [level]
  (when-not (and (number? level) (m/< 0.0 (double level) 1.0))
    (ic/invalid-input! :bootstrap :level (str "level has to be in (0, 1), got " level) {:value level :expected "(0, 1)"})))

(defrecord Bootstrap [fit fits ^long failed interval estimation]
  Fn ;; marker: `fn?` is true, as for the interpolators
  IFn
  (invoke [_ x] (fit x))
  (invoke [_ x level]
    (check-level! level)
    (let [t0 (double (fit x))]
      (if (empty? fits)
        [t0 t0]
        (let [[lo hi] (interval {:t0 t0 :ts (mapv #(double (% x)) fits)} (m/- 1.0 (double level)) estimation)]
          [(double lo) (double hi)])))))

(defn- prepared
  "`[xs ys]` of the pairs `[x y]`, sorted by x, with the mean y of repeated x."
  [pairs]
  (let [groups (reduce (fn [acc [x y]] (update acc (double x) (fnil conj []) (double y))) (sorted-map) pairs)]
    [(vec (keys groups))
     (mapv (fn [vs] (m// (double (reduce m/+ 0.0 vs)) (double (count vs)))) (vals groups))]))

(defn bootstrap
  "Bootstrap band of the fit `fit-fn` - `(fn [xs ys] f)` returning a function of x - of the data `xs`, `ys`.

  The data (and every bootstrap sample) is sorted by x, with the mean y of repeated x values, so methods needing
  strictly increasing nodes can be used. The `:samples` samples of the `(x, y)` pairs (with replacement) are drawn by
  [[fastmath.stats.bootstrap/bootstrap]]; a sample whose fit throws (for example too few distinct points for the
  method) is skipped and counted.

  Parameters:

  - `fit-fn` (function): `(fn [xs ys] f)`, `f` a function of one number
  - `xs`, `ys` (sequences of numbers): the data, of equal length
  - `opts` (map, optional):
      - `:samples` (positive integer, default: `200`): number of bootstrap samples
      - `:rng` (random number generator, default: a new unseeded generator, see [[fastmath.random/ensure-rng]]): seed it to reproduce the band
      - `:ci` (keyword, default: `:percentile`): the confidence interval of a band point, one of [[intervals]] (`:percentile`, `:basic`, `:normal`, `:bc`)
      - `:estimation` (keyword, default: `:legacy`): quantile estimation strategy of the interval (see [[fastmath.stats/quantiles]]; not used by `:normal`)

  Returns a `Bootstrap` record, a function: `(f x)` - the fit of the (prepared) data at `x`; `(f x level)` -
  `[lo hi]`, the `:ci` interval at confidence `level` in (0, 1) of the refitted values at `x` around the fit;
  `[fit fit]` when no sample could be refitted. Fields: `:fit`, `:fits` (the refitted functions), `:failed`
  (skipped samples). The fits are computed once; a band point costs `O(samples)`.

  Throws `ex-info` (`:type :invalid-input`, `:method :bootstrap`) for different lengths of `xs` and `ys`
  (`:length-mismatch`), a `:samples` that isn't a positive integer (`:samples`), an unknown `:ci` (`:ci`) or a
  `level` outside (0, 1) (`:level`); errors of the fit of the data itself propagate."
  ([fit-fn xs ys] (bootstrap fit-fn xs ys nil))
  ([fit-fn xs ys {:keys [samples rng ci estimation]}]
   (let [samples (if (nil? samples) default-samples samples)
         ci (or ci :percentile)
         n (count xs)]
     (when-not (m/== n (count ys))
       (ic/invalid-input! :bootstrap :length-mismatch
                          (str "xs and ys have different lengths: " n " and " (count ys))
                          {:xs-count n :ys-count (count ys)}))
     (when-not (and (integer? samples) (m/pos? (long samples)))
       (ic/invalid-input! :bootstrap :samples (str "samples has to be a positive integer, got " samples)
                          {:value samples :expected "a positive integer"}))
     (when-not (contains? intervals ci)
       (ic/invalid-input! :bootstrap :ci (str ":ci has to be one of " (vec (sort (keys intervals))) ", got " ci)
                          {:value ci :expected (set (keys intervals))}))
     (let [pairs (mapv vector xs ys)
           [pxs pys] (prepared pairs)
           fit (fit-fn pxs pys)
           drawn (:samples (boot/bootstrap {:data pairs} {:samples samples
                                                         :rng (r/ensure-rng rng)
                                                         :distribution :categorical-distribution}))
           refits (map (fn [sample] (let [[sxs sys] (prepared sample)]
                                      (try (fit-fn sxs sys) (catch Exception _ nil))))
                       drawn)
           fits (vec (remove nil? refits))]
       (->Bootstrap fit fits (m/- (long samples) (count fits)) (intervals ci) (or estimation :legacy))))))
