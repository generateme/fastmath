(ns fastmath.kernel.variogram
  "Semivariograms and empirical (experimental) estimators for kriging interpolation.

  Semivariogram models are created using `nugget`, `partial sill`, `range` and in some cases `beta` (exponent) parameters. All models can be fit to the `empirical` (experimental) semivariogram which can be created by using a variety of the estimators.

  Estimators calculate `gamma` semivariogram for given bin `h` (lag) (`n` - number of pairs, `diffs` - seq of differences)."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.distance :as dist]
            [fastmath.stats :as stats]
            [fastmath.optimization :as optim]
            [fastmath.random :as r])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn power
  "Power semivariogram model."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m/pow (m// x range) beta)
             (m/* psill)
             (m/+ nugget)))))

(defn exponential
  "Exponential semivariogram model"
  [{:keys [^double nugget ^double psill ^double range]}]
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m// x range)
             (m/-)
             (m/exp)
             (m/- 1.0)
             (m/* psill)
             (m/+ nugget)))))

(defn expower
  "Exponential-power semivariogram model."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m/pow (m// x range) beta)
             (m/-)
             (m/exp)
             (m/- 1.0)
             (m/* psill)
             (m/+ nugget)))))

(defn gaussian
  "Gaussian semivariogram model."
  [{:keys [^double nugget ^double psill ^double range]}]
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m// x range)
             (m/sq)
             (m/-)
             (m/exp)
             (m/- 1.0)
             (m/* psill)
             (m/+ nugget)))))

(defn spherical
  "Spherical semivariogram model."
  [{:keys [^double nugget ^double psill ^double range]}]
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (let [xr (m// x range)]
              (->> (m/cb xr)
                   (m/* 0.5)
                   (m/- (m/* 1.5 xr))
                   (m/* psill)
                   (m/+ nugget))))))

(defn cubic
  "Cubic semivariogram model."
  [{:keys [^double nugget ^double psill ^double range]}]
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (let [xr (m// x range)
                  xr2 (m/sq xr)
                  xr3 (m/* xr2 xr)
                  xr5 (m/* xr3 xr2)
                  xr7 (m/* xr5 xr2)]
              (->> (m/+ (m/* 7.0 xr2)
                        (m/* -8.75 xr3)
                        (m/* 3.5 xr5)
                        (m/* -0.75 xr7))
                   (m/* psill)
                   (m/+ nugget))))))

(defn pentaspherical
  "Pentaspherical semivariogram model."
  [{:keys [^double nugget ^double psill ^double range]}]
  (fn ^double [^double x]
    (cond
      (zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (let [xr (m// x range)
                  xr2 (m/* xr xr)]
              (->> (m/* xr2 0.375)
                   (m/+ -1.25)
                   (m/* xr2)
                   (m/+ 1.875)
                   (m/* psill xr)
                   (m/+ nugget))))))

(defn circular
  "Circular semivariogram model."
  [{:keys [^double nugget ^double psill ^double range]}]
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (let [xr (m// x range)]
              (->> (m/* m/TWO_INV_PI
                        (m/+ (m/* xr (m/sqrt (m/- 1.0 (m/* xr xr))))
                             (m/asin xr)))
                   (m/* psill)
                   (m/+ nugget))))))

(defn linear
  "Linear semivariogram model."
  [{:keys [^double nugget ^double psill ^double range]}]
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (->> (m// x range)
                 (m/* psill)
                 (m/+ nugget)))))

(defn hole
  "Hole effect semivariogram model."
  [{:keys [^double nugget ^double psill ^double range]}]
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m// x range)
             (m/sinc)
             (m/- 1.0)
             (m/* psill)
             (m/+ nugget)))))

(defn cauchy
  "Generalized Cauchy semivariogram model."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (-> (m// x range)
                 (m/sq)
                 (m/inc)
                 (m/pow beta)
                 (m//))
             (m/- 1.0)
             (m/* psill)
             (m/+ nugget)))))

(defn ->bessel
  "Creator of the Bessel (of the first kind) semivariogram model."
  [^double alpha]
  (let [g (m/gamma (m/inc alpha))]
    (fn [{:keys [^double nugget ^double psill ^double range ^double beta]}]
      (fn ^double [^double x]
        (if (m/zero? x) 0.0
            (->> (m/* g (m/pow (m// (m/* 2.0 range) x) beta)
                      (m/bessel-j alpha (m// x range)))
                 (m/- 1.0)
                 (m/* psill)
                 (m/+ nugget)))))))


(defn rbf->variogram
  "Build semivariogram model based on RBF kernel.

  Selected kernel should has a value `1.0` at zero and should be `0.0` at infinity. For example `:matern-c2` will work but `:thin-plate` will not. `range` parameter acts as reciprocal of RBF's `shape`."
  [kernel]
  (fn [{:keys [^double nugget ^double psill ^double range]}]
    (fn ^double [^double x]
      (if (m/zero? x) 0.0
          (->> (m// x range)
               ^double (kernel)
               (m/- 1.0)
               (m/* psill)
               (m/+ nugget))))))

(def ^:private semivariograms {:linear linear
                              :pentaspherical pentaspherical
                              :spherical spherical
                              :gaussian gaussian
                              :exponential exponential
                              :power power
                              :expower expower
                              :hole hole
                              :circular circular
                              :cubic cubic
                              :cauchy cauchy})

(def ^:private power-semivariograms #{:power :expower :cauchy :bessel})

;; Highly Robust Variogram Estimation
;; Marc G. Genton

;; Basic Steps in Geostatistics: The Variogram and Kriging
;; Margaret A. Oliver, Richard Webster

(defn matheron-estimator
  "Matheron (classical) empirical semivariogram estimator."
  ^double [_ diffs]
  (stats/mean (map m/sq diffs)))

(defn cressie-estimator
  "Cressie (robust) empirical semivariogram estimator."
  [^long n diffs]
  (-> (->> (map (fn [^double diff] (m/sqrt (m/abs diff))) diffs)
           (stats/mean))
      (m/fpow 4)
      (m// (m/+ 0.457 (m// 0.494 n) (m// 0.045 (m/* n n))))))

(defn highly-robust-estimator
  "Genton's higly robust empirical semivariogram estimator."
  [^long n diffs]
  (let [idiffs (map-indexed vector diffs)
        ds (for [[^long i1 ^double d1] idiffs
                 [^long i2 ^double d2] idiffs
                 :when (m/< i1 i2)]
             (m/abs (m/- d1 d2)))
        q (m// (m/combinations (inc (int (m/* 0.5 n))) 2)
               (count ds))]
    (m/sq (m/* 2.2191 (stats/quantile ds q)))))

(defn dowd-estimator
  "Dowd's empirical semivariogram estimator.

  It's the same as quantile estimator for `p` equal `0.5`."
  [_ diffs]
  (m/* 2.198109338317728 (stats/median (map m/sq diffs))))

(defn ->quantile-estimator
  "Create a quantile (Armstrong and Delfiner) empirical semivariogram estimator for given `p`."
  [^double p]
  (let [fact (m// (double (r/icdf (r/distribution :chi-squared) p)))]
    (fn [_ diffs]
      (m/* fact (stats/quantile (map m/sq diffs) p)))))

(defn- psi
  [^double t]
  (if (m/>= (m/abs t) 4.0)
    0.0
    (m/* t (m/sq (m/- 16.0 (m/sq t))))))

(defn- dpsi
  [^double t]
  (if (m/>= (m/abs t) 4.0)
    0.0
    (let [t2 (m/sq t)]
      (m/+ (m/* 5.0 (m/sq t2))
           (m/* -96.0 t2)
           256.0))))

;; ROBUST SEMIVARIOGRAM ESTIMATION 
;; IN THE PRESENCE OF INFLUENTIAL 
;; SPATIAL DATA VALUES
;; Richard F. Gunst and Molly I. Hartfield 

(defn robust-m-estimator
  "Robust M-estimator (Gunst and Hartfield) for empirical semivariogram."
  [^long n diffs]
  (let [adiffs (map (fn [^double d] (m/sqrt (m/abs d))) diffs)
        m (stats/median adiffs)
        s (m// (stats/median-absolute-deviation adiffs m) 0.6745)
        y (map (fn [^double az] (m// (m/- az m) s)) adiffs)
        spsi (stats/sum (map psi y))
        sdpsi (stats/sum (map dpsi y))]
    (m// (m/fpow (m/+ m (m/* s (m// spsi sdpsi))) 4)
         (m/+ 0.457 (m// 0.494 n) (m// 0.045 (m/* n n))))))

(defn bounding-box-diagonal
  "Length of the diagonal of bounding box of spatial points."
  ^double [xs]
  (let [mn (reduce v/emn xs)
        mx (reduce v/emx xs)]
    (v/dist mn mx)))

(defn remove-outliers-fence
  "Removes outliers using Tukey's fences (`k`=1.5)."
  [combined ys]
  (let [[^double q1 ^double q3] (stats/quantiles ys [0.25 0.75])
        iqr (m/* 1.5 (m/- q3 q1))
        lif-thr (m/- q1 iqr)
        uif-thr (m/+ q3 iqr)]
    (filter (fn [v] (m/<= lif-thr ^double (v 2) uif-thr)) combined)))

(defn remove-outliers-mad
  "Removes outliers using median absolute deviation."
  [combined ys]
  (let [m (stats/median ys)
        s (m// (stats/median-absolute-deviation ys m) 0.6745)]
    (filter (fn [v] (m/<= (m/abs (m// (m/- ^double (v 2) m) s)) 3.0)) combined)))

(defn- maybe-remove-outliers
  [combined ys remove-outliers?]
  (cond
    (= :mad remove-outliers?) (remove-outliers-mad combined ys)
    remove-outliers? (remove-outliers-fence combined ys)
    :else combined))

(defn empirical
  "Empirical (experimental) semivariogram.

  Arguments:
  * `xs` - positions
  * `ys` - values
  * parameters (optional):
      * `:bins` - number of bins, size of the semivariogram (default: `15`)
      * `:cutoff` - semivariogram cutoff (default: bounding box diagonal divided by `diagonal-den`)
      * `:diagonal-den` - denomiator of bounding box diagonal (default: `3`)
      * `:estimator` - estimator name or a function (default: `:classical`)
      * `:quantile` - quantile for quantile estimator (default: `0.9`)
      * `:remove-outliers?` - should outliers be removed from raw data? (default: `false`)

  Defined estimators are: `:classical`/`:matheron` (default), `:cressie`, `:genton`/`:highly-robust`, `:dowd`, `:quantile`, `:m-robust`.

  `:remove-outliers?` can be a `true` value for Tukey's fences criterion or `:mad` for median absolute deviation criterion.

  Function return a list of maps sorted by lag `h`, containing:

  * `:n` - number of points in given bin
  * `:h` - average lag
  * `:gamma` - semivariogram estimation"
  ([xs ys] (empirical xs ys nil))
  ([xs ys {:keys [^double cutoff ^double diagonal-den ^long bins estimator ^double quantile
                  remove-outliers?]
           :or {bins 15 diagonal-den 3.0 estimator :classical quantile 0.9
                remove-outliers? false}}]
   (let [cutoff (or cutoff (/ (bounding-box-diagonal xs) diagonal-den))
         combined (-> (map vector (range) xs ys)
                      (maybe-remove-outliers ys remove-outliers?))
         distances (->> (for [[^long id1 x1 ^double y1] combined
                              [^long id2 x2 ^double y2] combined
                              :when (< id1 id2)
                              :let [d (dist/euclidean x1 x2)]
                              :when (and (pos? d) (<= d cutoff))]
                          (Vec2. d (m/- y1 y2)))
                        (sort-by first))
         max-dist (.x ^Vec2 (last distances))
         splits (rest (m/slice-range 0 max-dist (inc bins)))
         estimator-fn (if (fn? estimator)
                        estimator
                        (case estimator
                          :classical matheron-estimator
                          :matheron matheron-estimator
                          :cressie cressie-estimator
                          :genton highly-robust-estimator
                          :highly-robust highly-robust-estimator
                          :dowd dowd-estimator
                          :quantile (->quantile-estimator quantile)
                          :m-robust robust-m-estimator))]
     (loop [buff []
            distances distances
            [^double c & rsplits] splits]
       (let [found (take-while (fn [^Vec2 v] (m/<= (.x v) c)) distances)]
         (if-not (seq found)
           (recur buff distances rsplits)
           (let [n (count found)
                 buff (conj buff {:n n
                                  :h (stats/mean (map first found))
                                  :gamma (m/* 0.5 ^double (estimator-fn n (map second found)))})]
             (if (seq rsplits)
               (recur buff (drop n distances) rsplits)
               buff))))))))

(defn- weights-n-v2  [ns vs] (map (fn [^long n ^double v] (m// n (m/max m/EPSILON (m/sq v)))) ns vs))
(defn- weights-n-v3  [ns vs] (map (fn [^long n ^double v] (m// n (m/max m/EPSILON (m/cb v)))) ns vs))
(defn- est-max ^double [^double a ^double b] (m/max a b))
(defn- est-sq  ^double [^double a ^double b] (m/sq (m/- a b)))
(defn- est-abs ^double [^double a ^double b] (m/abs (m/- a b)))

(defn- infer-target-args
  [defaults params]
  (->> params
       (reduce (fn [buff t]
                 (if-not (t defaults)
                   (conj buff t)
                   buff)) [])))

(defn fit-params
  "Fits a model and a pair of selected semivariogram model and fitted parameters."
  ([empirical-semivariogram semivariogram-model] (fit-params empirical-semivariogram semivariogram-model nil))
  ([empirical-semivariogram semivariogram-model {:keys [estimation weights defaults ^double order]
                                                 :or {estimation :sq order 0.0}}]
   (let [semivariogram-fn (cond
                            (fn? semivariogram-model) (rbf->variogram semivariogram-model)
                            (= :bessel semivariogram-model) (->bessel order)
                            :else (semivariograms semivariogram-model))
         
         target-args (infer-target-args defaults (if (power-semivariograms semivariogram-model)
                                                   [:nugget :psill :range :beta]
                                                   [:nugget :psill :range]))]

     (if (seq target-args)
       ;; fitting is needed
       (let [est-fn (case estimation :abs est-abs :max est-max est-sq)
             hs (map :h empirical-semivariogram)
             gammas (map :gamma empirical-semivariogram)
             ns (map :n empirical-semivariogram)

             max-h (last hs)         
             max-gamma (stats/maximum gammas)


             bounds-map {:nugget [0.0 max-gamma]
                         :psill [0.0 (m/* 10.0 max-gamma)]
                         :range [1.0 max-h]
                         :beta [-50.0 50.0]}
             
             bounds (map bounds-map target-args)

             target (condp = weights
                      :ngg (fn [& r] (let [sv (semivariogram-fn (merge (zipmap target-args r) defaults))
                                          svh (map sv hs)
                                          weights-seq (weights-n-v2 ns svh)]
                                      (-> (map (fn [^double gamma- ^double gamma]
                                                 (est-fn gamma- gamma)) svh gammas)
                                          (stats/wmean weights-seq))))

                      :ngg2 (let [ng (map m/fast* ns gammas)]
                              (fn [& r] (let [sv (semivariogram-fn (merge (zipmap target-args r) defaults))
                                             svh (map sv hs)
                                             weights-seq (weights-n-v3 ng svh)]
                                         (-> (map (fn [^double gamma- ^double gamma]
                                                    (est-fn gamma- gamma)) svh gammas)
                                             (stats/wmean weights-seq)))))
                      
                      (let [weights-seq (cond
                                          (not weights) (repeat (count empirical-semivariogram) 1.0)
                                          (sequential? weights) weights
                                          :else (case weights
                                                  :n ns
                                                  :nhh (weights-n-v2 ns hs)))]                    
                        (fn [& r] (let [sv (semivariogram-fn (merge (zipmap target-args r) defaults))]
                                   (-> (map (fn [^double h ^double gamma]
                                              (est-fn (sv h) gamma)) hs gammas)
                                       (stats/wmean weights-seq))))))
             
             m (optim/scan-and-minimize :lbfgsb target {:bounds bounds})]

         [semivariogram-fn (->> m
                                (first)
                                (zipmap target-args)
                                (merge defaults))])
       ;; no fitting needed
       [semivariogram-fn defaults]))))

(defn fit
  "Fits semivariogram model parameter to a empirical semivariogram data.

  Arguments:

  * `empirical-semivariogram` - empirical semivariogram data as returned by `empirical` function.
  * `semivariogram-model` - name of the semivariogram model
  * parameters:
      * `:estimation` - estimation used to fit the model, `:sq` - least squares (default), `:abs` - least absolute values
      * `:weights` - fitting weights (default: `nil`, no weights)
      * `:order` - order for `:bessel` semivariogram model
      * `:defaults` - a map containing semivariogram model parameters which should be a fixed value, fitting will be done for the lacking ones only.

  Weights can be a sequence of weights or one of the following methods:
  
  * `:n` - N, number of points in a bin
  * `:nhh` - N divided by a squared lag (h)
  * `:ngg` - N divided by a squared gamma (estimated from a model)
  * `:ngg2` -N divided by a cubed gamma (estimated from a model) and multiplied by empirical gamma.
  
  Fit is done by using L-BFGS-B numerical optimization.

  Returns fitted semivariogram model."
  ([empirical-variogram semivariogram-model] (fit empirical-variogram semivariogram-model nil))
  ([empirical-variogram semivariogram-model parameters]
   (let [[f p] (fit-params empirical-variogram semivariogram-model parameters)]
     (f p))))
