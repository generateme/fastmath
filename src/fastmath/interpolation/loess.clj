(ns fastmath.interpolation.loess
  "LOESS (locally weighted linear regression) smoothing followed by cubic spline interpolation of the smoothed values."
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic]
            [fastmath.interpolation.cubic :as cubic])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^:private ^:const default-bandwidth 0.3)
(def ^:private ^:const default-iters 2)
(def ^:private ^:const default-accuracy 1.0e-12)

(defn- window-size
  "Number of nearest neighbours used in every local fit: the integer part of `bandwidth * n`, rounded up when
  the product is within 1e-10 below an integer (statsmodels rule)."
  ^long [^double bandwidth ^long n]
  (long (m/+ (m/* bandwidth (double n)) 1.0e-10)))

(defn- tricube
  "Tricube weight of the relative distance `d` in [0, 1]."
  ^double [^double d]
  (let [d3 (m/* d d d)
        u (m/- 1.0 d3)]
    (m/* u u u)))

(defn- fit-pass
  "One pass of local linear regressions at every node. `w` are the weights of the observations (user weights
  times robustness weights), `k` is the number of the nearest neighbours.

  The window of `k` consecutive nodes slides to the right while the evaluated node is closer to the next
  node than to the first node of the window, so it always contains the `k` nearest neighbours. Weights are
  tricube functions of the distance divided by the largest distance in the window. When fewer than two
  weights are positive the observation itself is kept."
  [^doubles xs ^doubles ys ^doubles w k]
  (let [n (alength xs)
        k (long k)
        yfit (double-array n)
        wt (double-array n)]
    (loop [i (long 0)
           left (long 0)]
      (when (m/< i n)
        (let [xv (Array/aget xs i)
              left (long (loop [l left]
                           (if (and (m/< (m/+ l k) n)
                                    (m/> xv (m/* 0.5 (m/+ (Array/aget xs l) (Array/aget xs (m/+ l k))))))
                             (recur (m/inc l))
                             l)))
              right (m/+ left k)
              radius (m/max (m/- xv (Array/aget xs left)) (m/- (Array/aget xs (m/dec right)) xv))]
          (loop [j left
                 sw 0.0
                 cnt (long 0)]
            (if (m/< j right)
              (let [wj (m/* (tricube (m// (m/abs (m/- (Array/aget xs j) xv)) radius))
                            (Array/aget w j))]
                (Array/aset wt j wj)
                (recur (m/inc j) (m/+ sw wj) (if (m/> wj 1.0e-12) (m/inc cnt) cnt)))
              (if (m/< cnt 2)
                (Array/aset yfit i (Array/aget ys i))
                (let [swx (double (loop [j left
                                         s 0.0]
                                    (if (m/< j right)
                                      (recur (m/inc j) (m/+ s (m/* (m// (Array/aget wt j) sw) (Array/aget xs j))))
                                      s)))
                      wsq (m/max (double (loop [j left
                                                s 0.0]
                                           (if (m/< j right)
                                             (let [dx (m/- (Array/aget xs j) swx)]
                                               (recur (m/inc j) (m/+ s (m/* (m// (Array/aget wt j) sw) dx dx))))
                                             s)))
                                 (m/* 1.0e-12 radius radius))
                      dxv (m/- xv swx)]
                  (Array/aset yfit i
                              (double (loop [j left
                                             s 0.0]
                                        (if (m/< j right)
                                          (let [wj (m// (Array/aget wt j) sw)
                                                pj (m/* wj (m/+ 1.0 (m// (m/* dxv (m/- (Array/aget xs j) swx)) wsq)))]
                                            (recur (m/inc j) (m/+ s (m/* pj (Array/aget ys j)))))
                                          s))))))))
          (recur (m/inc i) left))))
    yfit))

(defn- median-abs-residual
  "True median of absolute residuals (mean of the two central values for an even count)."
  ^double [^doubles ys ^doubles yfit]
  (let [n (alength ys)
        r (double-array n)]
    (dotimes [i n]
      (Array/aset r i (m/abs (m/- (Array/aget ys i) (Array/aget yfit i)))))
    (java.util.Arrays/sort r)
    (if (m/even? n)
      (m/* 0.5 (m/+ (Array/aget r (m/dec (m/quot n 2))) (Array/aget r (m/quot n 2))))
      (Array/aget r (m/quot n 2)))))

(defn- robustness-weights
  "Bisquare robustness weights of residuals scaled by 6 times their median. When the median is zero only
  exactly zero residuals keep the full weight."
  ^doubles [^doubles ys ^doubles yfit med]
  (let [n (alength ys)
        med (double med)
        w (double-array n)]
    (dotimes [i n]
      (let [r (m/abs (m/- (Array/aget ys i) (Array/aget yfit i)))
            u (if (m/zero? med)
                (if (m/pos? r) 1.0 0.0)
                (m/min 1.0 (m// r (m/* 6.0 med))))
            v (m/- 1.0 (m/* u u))]
        (Array/aset w i (m/* v v))))
    w))

(defn- smooth
  "Smoothed values at the nodes: local linear regression with `k` nearest neighbours, followed by `iters`
  robustness iterations. The iterations stop earlier when the median absolute residual is below `accuracy`."
  ^doubles [^doubles xs ^doubles ys ^doubles uw opts]
  (let [{:keys [k iters accuracy]} opts
        k (long k)
        iters (long iters)
        accuracy (double accuracy)
        n (alength xs)
        w (double-array n)]
    (loop [it (long 0)
           rw nil]
      (dotimes [i n]
        (Array/aset w i (m/* (if uw (Array/aget uw i) 1.0)
                             (if rw (Array/aget ^doubles rw i) 1.0))))
      (let [yfit (fit-pass xs ys w k)]
        (if (m/< it iters)
          (let [med (median-abs-residual ys yfit)]
            (if (m/< med accuracy)
              yfit
              (recur (m/inc it) (robustness-weights ys yfit med))))
          yfit)))))

(defn- check-options!
  [n bandwidth iters accuracy weights]
  (when-not (and (number? bandwidth) (m/> (double bandwidth) 0.0) (m/<= (double bandwidth) 1.0))
    (ic/invalid-input! :loess :bandwidth
                       (str "bandwidth has to be in (0, 1], got " bandwidth)
                       {:value bandwidth :expected "(0, 1]"}))
  (let [k (window-size (double bandwidth) (long n))]
    (when (m/< k 2)
      (ic/invalid-input! :loess :bandwidth
                         (str "bandwidth " bandwidth " with " n " points gives a window of " k
                              " nearest neighbours, at least 2 are required (bandwidth * n >= 2)")
                         {:value bandwidth :count n :window k :min-window 2})))
  (when-not (and (number? iters) (m/>= (double iters) 0.0))
    (ic/invalid-input! :loess :iters (str "iters has to be a non-negative number, got " iters)
                       {:value iters :expected ">= 0"}))
  (when-not (and (number? accuracy) (m/>= (double accuracy) 0.0))
    (ic/invalid-input! :loess :accuracy (str "accuracy has to be a non-negative number, got " accuracy)
                       {:value accuracy :expected ">= 0"}))
  (when weights
    (let [^doubles w weights]
      (when-not (m/== (alength w) (long n))
        (ic/invalid-input! :loess :length-mismatch
                           (str "weights and xs have different lengths: " (alength w) " and " n)
                           {:weights-count (alength w) :xs-count n}))
      (dotimes [i (alength w)]
        (let [wi (Array/aget w i)]
          (when-not (and (m/valid-double? wi) (m/>= wi 0.0))
            (ic/invalid-input! :loess :invalid-weights
                               (str "weights have to be finite and non-negative, got " wi " at index " i)
                               {:index i :value wi})))))))

(defn loess
  "Creates a one dimensional LOESS (locally weighted linear regression) smoother returning an interpolating function.

  Every node is smoothed by a weighted linear regression over its `k` nearest neighbours with tricube weights, optionally repeated with robustness (bisquare) weights that downweight outliers. A natural cubic spline is then fitted through the smoothed values, so the returned function is continuous and smooth and can be evaluated anywhere. The smoothed values equal `statsmodels.nonparametric.lowess(frac=bandwidth, it=iters, delta=0)`.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing, at least 2 points
  - `ys` (sequence of numbers): values at the nodes, the same count as `xs`
  - `opts` (map, optional):
      - `:bandwidth` (number in (0, 1], default: `0.3`): fraction of nodes used in each local fit. The window holds `k = floor(bandwidth * n + 1e-10)` true nearest neighbours of the evaluated node (the node itself included). `k` has to be at least 2, which means `bandwidth * n >= 2`; for the default `0.3` at least 7 points are needed, for `1.0` at least 2 points. `k = 2` returns the data unchanged, because the farthest neighbour gets the zero weight.
      - `:iters` (non-negative integer, default: `2`): number of robustness iterations. `0` is the plain LOESS. In every iteration the residuals are divided by 6 times their true median (the mean of the two central values for an even count) and the bisquare function gives the weights; points with residual of at least 6 medians have zero weight.
      - `:accuracy` (non-negative number, default: `1.0e-12`): the iterations stop when the median of absolute residuals falls below this value (for example when the data lie exactly on a line). With `0.0` the iterations always run, as statsmodels does.
      - `:weights` (sequence of non-negative numbers, default: `nil`): weights of the observations, the same count as `xs`. They are multiplied into the tricube weights (and into the robustness weights) of every local fit; an observation with the weight zero is ignored. `nil` means equal weights.

  Returns a function `f` of one argument `x` returning the interpolated `double` value; it is the natural cubic spline through the smoothed values (beyond the domain the end cubic polynomials are continued). The curve does not pass through the data points unless `k` is 2 or the data is exactly linear.

  A local fit with less than two positive weights keeps the original value. The variance of the weighted regression abscissas is bounded from below by 1e-12 times the squared window radius, so the result does not depend on the scale of `xs`. Each pass costs `O(n k)`.

  Throws `ex-info` (with `:type` `:invalid-input`, `:method` `:loess` and a `:property`) when the lengths differ (`:length-mismatch`), there are fewer than 2 points (`:too-few-points`), `xs` is not finite (`:not-finite`) or not strictly increasing (`:not-increasing`), `:bandwidth` is outside (0, 1] or gives `k < 2` (`:bandwidth`), `:iters` or `:accuracy` is negative (`:iters`, `:accuracy`) or `:weights` have a wrong length (`:length-mismatch`), a negative or a non-finite value (`:invalid-weights`).

  See also [[cubic-smoothing]], [[cubic]], [[interpolation]]."
  ([xs ys] (loess xs ys nil))
  ([xs ys {:keys [bandwidth iters accuracy weights]}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         uw (when weights (m/seq->double-array weights))
         bandwidth (if (nil? bandwidth) default-bandwidth bandwidth)
         iters (if (nil? iters) default-iters iters)
         accuracy (if (nil? accuracy) default-accuracy accuracy)]
     (ic/validate-xy! :loess xs ys {:min-points 2})
     (check-options! (alength xs) bandwidth iters accuracy uw)
     (let [opts {:k (window-size (double bandwidth) (alength xs))
                 :iters (long iters)
                 :accuracy (double accuracy)}]
       (cubic/cubic xs (smooth xs ys uw opts))))))
