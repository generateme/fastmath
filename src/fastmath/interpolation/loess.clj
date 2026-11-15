(ns fastmath.interpolation.loess
  "LOESS (locally weighted linear or quadratic regression) smoothing followed by cubic spline interpolation of the
  smoothed values."
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic]
            [fastmath.interpolation.cubic :as cubic])
  (:import [fastmath.java Array]
           [clojure.lang Fn IFn]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^:private ^:const default-bandwidth 0.3)
(def ^:private ^:const default-degree 1)
(def ^:private ^:const default-iters 2)
(def ^:private ^:const default-accuracy 1.0e-12)

(defn- window-size
  "Number of nearest neighbours used in every local fit: the integer part of `bandwidth * n`, rounded up when
  the product is within 1e-10 below an integer (statsmodels rule)."
  ^long [^double bandwidth ^long n]
  (long (m/+ (m/* bandwidth n) 1.0e-10)))

(defn- tricube
  "Tricube weight of the relative distance `d` in [0, 1]."
  ^double [^double d]
  (let [d3 (m/* d d d)
        u (m/- 1.0 d3)]
    (m/* u u u)))

(defrecord ^:no-doc LocalFits [^doubles fit ^doubles hat-diag ^doubles hat-sq ^doubles hat-var])

(defn- linear-row!
  "Row of the hat matrix of the local linear fit at `xv` over the window `[left, right)` into `p` (`p[j - left]`):
  `(wt_j / sw) * (1 + (xv - mean) (x_j - mean) / var)`, `mean` and `var` the weighted mean and variance of the
  window (the variance bounded below by 1e-12 times the squared `radius`)."
  [^doubles p ^doubles xs ^doubles wt sw left right xv radius]
  (let [sw (double sw) left (long left) right (long right) xv (double xv) radius (double radius)
        swx (double (loop [j left s 0.0]
                      (if (m/< j right)
                        (recur (m/inc j) (m/+ s (m/* (m// (Array/aget wt j) sw) (Array/aget xs j))))
                        s)))
        wsq (m/max (double (loop [j left s 0.0]
                             (if (m/< j right)
                               (let [dx (m/- (Array/aget xs j) swx)]
                                 (recur (m/inc j) (m/+ s (m/* (m// (Array/aget wt j) sw) dx dx))))
                               s)))
                   (m/* 1.0e-12 radius radius))
        dxv (m/- xv swx)]
    (loop [j left]
      (when (m/< j right)
        (Array/aset p (m/- j left) (m/* (m// (Array/aget wt j) sw)
                                        (m/+ 1.0 (m// (m/* dxv (m/- (Array/aget xs j) swx)) wsq))))
        (recur (m/inc j))))))

(def ^{:private true :const true :tag 'double} quadratic-det-tolerance 1.0e-12)

(defn- quadratic-row!
  "Row of the hat matrix of the local quadratic fit at `xv` over the window `[left, right)` into `p`
  (`p[j - left]`): `(wt_j / sw) * (c0 + c1 u_j + c2 u_j^2)`, `u_j = (x_j - xv) / radius`, `c` - the first row of the
  inverse of the weighted moment matrix of `1, u, u^2`. Returns false, `p` untouched, when the moment matrix is
  (nearly) singular - fewer than three distinct points carry weight."
  [^doubles p ^doubles xs ^doubles wt sw left right xv radius]
  (let [sw (double sw) left (long left) right (long right) xv (double xv) radius (double radius)
        moments (double-array 5)]
    (loop [j left]
      (when (m/< j right)
        (let [wj (m// (Array/aget wt j) sw)
              u (m// (m/- (Array/aget xs j) xv) radius)
              u2 (m/* u u)]
          (Array/aset moments 0 (m/+ (Array/aget moments 0) wj))
          (Array/aset moments 1 (m/+ (Array/aget moments 1) (m/* wj u)))
          (Array/aset moments 2 (m/+ (Array/aget moments 2) (m/* wj u2)))
          (Array/aset moments 3 (m/+ (Array/aget moments 3) (m/* wj u2 u)))
          (Array/aset moments 4 (m/+ (Array/aget moments 4) (m/* wj u2 u2))))
        (recur (m/inc j))))
    (let [m0 (Array/aget moments 0) m1 (Array/aget moments 1) m2 (Array/aget moments 2)
          m3 (Array/aget moments 3) m4 (Array/aget moments 4)
          c00 (m/- (m/* m2 m4) (m/* m3 m3))
          c01 (m/- (m/* m2 m3) (m/* m1 m4))
          c02 (m/- (m/* m1 m3) (m/* m2 m2))
          det (m/+ (m/* m0 c00) (m/+ (m/* m1 c01) (m/* m2 c02)))]
      (if (m/<= det quadratic-det-tolerance)
        false
        (let [c0 (m// c00 det) c1 (m// c01 det) c2 (m// c02 det)]
          (loop [j left]
            (when (m/< j right)
              (let [u (m// (m/- (Array/aget xs j) xv) radius)]
                (Array/aset p (m/- j left) (m/* (m// (Array/aget wt j) sw) (m/+ c0 (m/* u (m/+ c1 (m/* c2 u)))))))
              (recur (m/inc j))))
          true)))))

(defn- fit-pass
  "One pass of local regressions (linear, or quadratic for `degree` 2) at every node. `w` are the weights of the
  observations (user weights times robustness weights), `k` is the number of the nearest neighbours.

  The window of `k` consecutive nodes slides to the right while the evaluated node is closer to the next
  node than to the first node of the window, so it always contains the `k` nearest neighbours. Weights are
  tricube functions of the distance divided by the largest distance in the window. When fewer than two
  weights are positive the observation itself is kept; a quadratic fit with fewer than three distinct weighted
  points is linear.

  The fit at node `i` is `sum_j L[i,j] * y_j`, `L` - the hat matrix of the pass (rows from [[linear-row!]] /
  [[quadratic-row!]]). Returns `LocalFits`: `:fit` - the fitted values and, with `hat?`, per node `i`, `:hat-diag` -
  `L[i,i]`, `:hat-sq` - `sum_j L[i,j]^2` and `:hat-var` - `sum_j L[i,j]^2 / uw_j` (`uw` - the user weights, `nil` for
  equal weights; a zero weight has zero `L[i,j]`); without `hat?` these are `nil`."
  [^doubles xs ^doubles ys ^doubles w ^doubles uw k degree hat?]
  (let [n (alength xs)
        k (long k)
        quadratic? (m/== 2 (long degree))
        yfit (double-array n)
        wt (double-array n)
        row (double-array k)
        l-diag (when hat? (double-array n))
        l-sq (when hat? (double-array n))
        l-var (when hat? (double-array n))]
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
                (do (Array/aset yfit i (Array/aget ys i))
                    (when hat?
                      (Array/aset l-diag i 1.0)
                      (Array/aset l-sq i 1.0)
                      (Array/aset l-var i (if uw (m// 1.0 (Array/aget uw i)) 1.0))))
                (do
                  (when-not (and quadratic? (m/>= cnt 3) (quadratic-row! row xs wt sw left right xv radius))
                    (linear-row! row xs wt sw left right xv radius))
                  (Array/aset yfit i
                              (double (loop [j left
                                             s 0.0]
                                        (if (m/< j right)
                                          (recur (m/inc j) (m/+ s (m/* (Array/aget row (m/- j left)) (Array/aget ys j))))
                                          s))))
                  (when hat?
                    (loop [j left
                           sq 0.0
                           vr 0.0]
                      (if (m/< j right)
                        (let [pj (Array/aget row (m/- j left))
                              pj2 (m/* pj pj)]
                          (when (m/== j i) (Array/aset l-diag i pj))
                          (recur (m/inc j) (m/+ sq pj2)
                                 (if (m/zero? pj2) vr (m/+ vr (if uw (m// pj2 (Array/aget uw j)) pj2)))))
                        (do (Array/aset l-sq i sq)
                            (Array/aset l-var i vr)))))))))
          (recur (m/inc i) left))))
    (->LocalFits yfit l-diag l-sq l-var)))

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
  "Smoothed values at the nodes: local regression of `degree` (1 or 2) with `k` nearest neighbours, followed by `iters`
  robustness iterations. The iterations stop earlier when the median absolute residual is below `accuracy`.
  Returns the `LocalFits` of the final pass (see `fit-pass`)."
  [^doubles xs ^doubles ys ^doubles uw opts]
  (let [{:keys [k degree iters accuracy]} opts
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
      ;; the hat matrix rows only in the final pass: the last iteration, or a pass the accuracy stops (redone)
      (let [last? (m/>= it iters)
            fits (fit-pass xs ys w uw k degree last?)
            ^doubles yfit (:fit fits)]
        (if last?
          fits
          (let [med (median-abs-residual ys yfit)]
            (if (m/< med accuracy)
              (fit-pass xs ys w uw k degree true)
              (recur (m/inc it) (robustness-weights ys yfit med)))))))))

(defrecord Loess [spline stderr-spline ^double sigma ^double df ^double trace]
  Fn ;; marker: `fn?` stays true, as for the plain function loess returned before
  IFn
  (invoke [_ x] (spline x))
  (invoke [_ x stderr?]
    (if stderr?
      [(double (spline x)) (double (stderr-spline x))]
      (spline x))))

(def ^{:private true :const true :tag 'double} df-tolerance 1.0e-10)

(defn- model
  "The `Loess` record of the final pass `fits` (`LocalFits`, see `fit-pass`) at the nodes `xs` (data `ys`, user
  weights `uw` or `nil`): residual scale from `delta1 = n - 2 tr(L) + tr(L^T L)` (observations with the weight zero
  left out), standard errors at the nodes `sigma * sqrt(sum_j L[i,j]^2 / uw_j)`, a natural cubic spline through
  both (standard errors `NaN` outside the nodes' range)."
  [^doubles xs ^doubles ys ^doubles uw fits]
  (let [^doubles yfit (:fit fits)
        ^doubles l-diag (:hat-diag fits)
        ^doubles l-sq (:hat-sq fits)
        ^doubles l-var (:hat-var fits)
        n (alength xs)
        sums (double-array 4) ;; included count, tr(L), tr(L^T L), weighted RSS
        _ (dotimes [i n]
            (when (or (nil? uw) (m/pos? (Array/aget uw i)))
              (let [r (m/- (Array/aget ys i) (Array/aget yfit i))]
                (Array/aset sums 0 (m/inc (Array/aget sums 0)))
                (Array/aset sums 1 (m/+ (Array/aget sums 1) (Array/aget l-diag i)))
                (Array/aset sums 2 (m/+ (Array/aget sums 2) (Array/aget l-sq i)))
                (Array/aset sums 3 (m/+ (Array/aget sums 3) (m/* (if uw (Array/aget uw i) 1.0) r r))))))
        trace (Array/aget sums 1)
        rss (Array/aget sums 3)
        delta1 (m/+ (m/- (Array/aget sums 0) (m/* 2.0 trace)) (Array/aget sums 2))
        df (if (m/> delta1 df-tolerance) delta1 0.0)
        sigma (if (m/pos? df) (m/sqrt (m// rss df)) ##NaN)
        spline (cubic/cubic xs yfit)
        lo (Array/aget xs 0)
        hi (Array/aget xs (m/dec n))
        stderr-spline (if (m/nan? sigma)
                        (constantly ##NaN)
                        (let [c (cubic/cubic xs (map (fn [^long i] (m/* sigma (m/sqrt (Array/aget l-var i)))) (range n)))]
                          (fn ^double [x]
                            (let [x (double x)]
                              (if (and (m/<= lo x) (m/<= x hi)) (m/max 0.0 (double (c x))) ##NaN)))))]
    (->Loess spline stderr-spline sigma df trace)))

(defn- check-options!
  [n bandwidth degree iters accuracy weights]
  (when-not (contains? #{1 2} degree)
    (ic/invalid-input! :loess :degree (str "degree has to be 1 or 2, got " degree) {:value degree :expected #{1 2}}))
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
  "Creates a one dimensional LOESS (locally weighted linear or quadratic regression) smoother returning an interpolating function.

  Every node is smoothed by a weighted linear (or quadratic, `:degree 2`) regression over its `k` nearest neighbours with tricube weights, optionally repeated with robustness (bisquare) weights that downweight outliers. A natural cubic spline is then fitted through the smoothed values, so the returned function is continuous and smooth and can be evaluated anywhere. With `:degree 1` the smoothed values equal `statsmodels.nonparametric.lowess(frac=bandwidth, it=iters, delta=0)`; with `:degree 2` and `:iters 0` R `loess(span = bandwidth, degree = 2, family = \"gaussian\", surface = \"direct\")` at the nodes.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing, at least 2 points
  - `ys` (sequence of numbers): values at the nodes, the same count as `xs`
  - `opts` (map, optional):
      - `:bandwidth` (number in (0, 1], default: `0.3`): fraction of nodes used in each local fit. The window holds `k = floor(bandwidth * n + 1e-10)` true nearest neighbours of the evaluated node (the node itself included). `k` has to be at least 2, which means `bandwidth * n >= 2`; for the default `0.3` at least 7 points are needed, for `1.0` at least 2 points. `k = 2` returns the data unchanged, because the farthest neighbour gets the zero weight.
      - `:degree` (`1` or `2`, default: `1`): degree of the local polynomials - local linear (statsmodels `lowess`) or local quadratic (R `loess` and ggplot2 `geom_smooth` default, `degree = 2`). A quadratic fit whose window has fewer than three distinct weighted points is linear.
      - `:iters` (non-negative integer, default: `2`): number of robustness iterations. `0` is the plain LOESS. In every iteration the residuals are divided by 6 times their true median (the mean of the two central values for an even count) and the bisquare function gives the weights; points with residual of at least 6 medians have zero weight.
      - `:accuracy` (non-negative number, default: `1.0e-12`): the iterations stop when the median of absolute residuals falls below this value (for example when the data lie exactly on a line). With `0.0` the iterations always run, as statsmodels does.
      - `:weights` (sequence of non-negative numbers, default: `nil`): weights of the observations, the same count as `xs`. They are multiplied into the tricube weights (and into the robustness weights) of every local fit; an observation with the weight zero is ignored. `nil` means equal weights.

  Returns a `Loess` record, a function: `(f x)` returns the interpolated `double` value - the natural cubic spline through the smoothed values (beyond the domain the end cubic polynomials are continued); the curve does not pass through the data points unless `k` is 2 or the data is exactly linear. `(f x true)` returns `[fit stderr]`, the standard error of the fit (as R `predict.loess(se = TRUE)`, with the hat matrix `L` of the final pass, robustness weights included): `sigma * sqrt(sum_j L[i,j]^2 / w_j)` at the nodes (`w` - the user weights, 1 without them) and a natural cubic spline through these values between them (clamped at 0); outside the range of the nodes the standard error is `NaN` (no information there, as R). The record holds `:sigma` - the residual scale `sqrt(sum w_i r_i^2 / delta1)`, `:df` - `delta1 = n - 2 tr(L) + tr(L^T L)` (R's `one.delta`, the degrees of freedom of a t interval) and `:trace` - `tr(L)`; observations with the weight zero are left out of the sums. With robustness iterations (`:iters` > 0) the hat matrix includes the robustness weights but `:sigma` uses the plain residuals, so outliers the robust fit downweights still widen the standard errors: the bands are conservative. With `k = 2` (or no residual degrees of freedom) `:df` is 0 and the standard errors are `NaN`; exactly linear data gives the standard error 0.

  A local fit with less than two positive weights keeps the original value. The variance of the weighted regression abscissas is bounded from below by 1e-12 times the squared window radius, so the result does not depend on the scale of `xs`. Each pass costs `O(n k)`.

  Throws `ex-info` (with `:type` `:invalid-input`, `:method` `:loess` and a `:property`) when the lengths differ (`:length-mismatch`), there are fewer than 2 points (`:too-few-points`), `xs` is not finite (`:not-finite`) or not strictly increasing (`:not-increasing`), `:bandwidth` is outside (0, 1] or gives `k < 2` (`:bandwidth`), `:degree` isn't 1 or 2 (`:degree`), `:iters` or `:accuracy` is negative (`:iters`, `:accuracy`) or `:weights` have a wrong length (`:length-mismatch`), a negative or a non-finite value (`:invalid-weights`).

  See also [[cubic-smoothing]], [[cubic]], [[interpolation]]."
  ([xs ys] (loess xs ys nil))
  ([xs ys {:keys [bandwidth degree iters accuracy weights]}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         uw (when weights (m/seq->double-array weights))
         bandwidth (if (nil? bandwidth) default-bandwidth bandwidth)
         degree (if (nil? degree) default-degree degree)
         iters (if (nil? iters) default-iters iters)
         accuracy (if (nil? accuracy) default-accuracy accuracy)]
     (ic/validate-xy! :loess xs ys {:min-points 2})
     (check-options! (alength xs) bandwidth degree iters accuracy uw)
     (let [opts {:k (window-size (double bandwidth) (alength xs))
                 :degree (long degree)
                 :iters (long iters)
                 :accuracy (double accuracy)}]
       (model xs ys uw (smooth xs ys uw opts))))))
