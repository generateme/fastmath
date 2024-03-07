(ns fastmath.kernel.variogram
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.distance :as dist]
            [fastmath.stats :as stats]
            [fastmath.optimization :as optim])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)

(defn power
  [{:keys [^double nugget ^double sill ^double range ^double beta]}]
  (fn ^double [^double x]
    (m/+ nugget (m/* sill (m/pow (m// x range) beta)))))

(defn exponential
  [{:keys [^double nugget ^double sill ^double range]}]
  (fn ^double [^double x]
    (m/+ nugget (m/* sill (m/- 1.0 (m/exp (m/- (m// x range))))))))

(defn expower
  [{:keys [^double nugget ^double sill ^double range ^double beta]}]
  (fn ^double [^double x]
    (m/+ nugget (m/* sill (m/- 1.0 (m/exp (m/- (m/pow (m// x range) beta))))))))

(defn gaussian
  [{:keys [^double nugget ^double sill ^double range]}]
  (fn ^double [^double x]
    (m/+ nugget (m/* sill (m/- 1.0 (m/exp (m/- (m/sq (m// x range)))))))))

(defn spherical
  [{:keys [^double nugget ^double sill ^double range]}]
  (fn ^double [^double x]
    (if (m/< range x)
      (m/+ nugget sill)
      (let [xr (m// x range)]
        (m/+ nugget (m/* sill (m/- (m/* 1.5 xr) (m/* 0.5 (m/cb xr)))))))))

(defn pentaspherical
  [{:keys [^double nugget ^double sill ^double range]}]
  (fn ^double [^double x]
    (cond
      (zero? x) nugget
      (m/< range x) (m/+ nugget sill)
      :else (let [xr (m// x range)
                  xr2 (m/* xr xr)]
              (->> (m/* xr2 0.375)
                   (m/+ -1.25)
                   (m/* xr2)
                   (m/+ 1.875)
                   (m/* sill xr)
                   (m/+ nugget))))))

(defn circular
  [{:keys [^double nugget ^double sill ^double range]}]
  (fn ^double [^double x]
    (cond
      (zero? x) nugget
      (m/< range x) (m/+ nugget sill)
      :else (let [xr (m// x range)]
              (->> (m/* m/TWO_INV_PI
                        (m/+ (m/* xr (m/sqrt (m/- 1.0 (m/* xr xr))))
                             (m/asin xr)))
                   (m/* sill)
                   (m/+ nugget))))))

(defn linear
  [{:keys [^double nugget ^double sill ^double range]}]
  (fn ^double [^double x]
    (if (m/< range x)
      (m/+ nugget sill)
      (let [xr (m// x range)]
        (m/+ nugget (m/* sill xr))))))

(defn hole
  [{:keys [^double nugget ^double sill ^double range]}]
  (fn ^double [^double x]
    (m/+ nugget (m/* sill (m/- 1.0 (m/sinc (m// x range)))))))

(defn rbf->variogram
  [kernel]
  (fn [{:keys [^double nugget ^double sill ^double range]}]
    (fn ^double [^double x]
      (m/+ nugget (m/* sill (m/- 1.0 ^double (kernel (m// x range))))))))

(def semi-variograms {:linear linear
                    :pentaspherical pentaspherical
                    :spherical spherical
                    :gaussian gaussian
                    :exponential exponential
                    :power power
                    :expower expower
                    :hole hole
                    :circular circular})

;; Highly Robust Variogram Estimation
;; Marc G. Genton

(defn matheron-estimator
  ^double [_ diffs]
  (stats/mean (map (comp m/sq) diffs)))

(defn cressie-estimator
  [^long n diffs]
  (-> (->> (map (fn [^double diff] (m/sqrt (m/abs diff))) diffs)
           (stats/mean))
      (m/fpow 4)
      (m// (m/+ 0.457 (m// 0.494 n)))))

(defn highly-robust-estimator
  [^long n diffs]
  (let [idiffs (map-indexed vector diffs)
        ds (for [[^long i1 ^double d1] idiffs
                 [^long i2 ^double d2] idiffs
                 :when (m/< i1 i2)]
             (m/abs (m/- d1 d2)))
        q (m// (m/combinations (inc (int (m/* 0.5 n))) 2)
               (count ds))]
    (m/sq (m/* 2.2191 (stats/quantile ds q)))))

(defn bounding-box-diagonal
  ^double [xs]
  (let [mn (reduce v/emn xs)
        mx (reduce v/emx xs)]
    (v/dist mn mx)))

(defn empirical
  ([xs ys] (empirical xs ys nil))
  ([xs ys {:keys [^double cutoff ^double diagonal-den ^long size estimator]
           :or {size 15 diagonal-den 3.0 estimator :classical}}]
   (let [cutoff (or cutoff (/ (bounding-box-diagonal xs) diagonal-den))
         combined (map vector (range) xs ys)
         distances (->> (for [[^long id1 x1 ^double y1] combined
                              [^long id2 x2 ^double y2] combined
                              :when (< id1 id2)
                              :let [d (dist/euclidean x1 x2)]
                              :when (and (pos? d) (<= d cutoff))]
                          (Vec2. d (m/- y1 y2)))
                        (sort-by first))
         max-dist (.x ^Vec2 (last distances))
         splits (rest (m/slice-range 0 max-dist (inc size)))
         estimator-fn (case estimator
                        :classical matheron-estimator
                        :matheron matheron-estimator
                        :cressie cressie-estimator
                        :highly-robust highly-robust-estimator)]
     (loop [buff []
            distances distances
            [^double c & rsplits] splits]
       (let [found (take-while (fn [^Vec2 v] (m/<= (.x v) c)) distances)]
         (if-not (seq found)
           (recur buff distances rsplits)
           (let [n (count found)
                 buff (conj buff {:n n
                                  :h (stats/mean (map first found))
                                  :gamma (m/* 0.5 (estimator-fn n (map second found)))})]
             (if (seq rsplits)
               (recur buff (drop n distances) rsplits)
               buff))))))))

(defn objective-least-squares
  ^double [variogram empirical-variogram]
  (->> empirical-variogram
       (map (fn [{:keys [^double gamma ^double h]}]
              (m/sq (m/- gamma (variogram h)))))
       (reduce m/fast+)))

(defn fit-params
  [empirical-variogram semi-variogram]
  (let [max-h (:h (last empirical-variogram))
        max-gamma (stats/maximum (map :gamma empirical-variogram))
        bounds [[0.0 max-gamma] [0.0 (m/* 10.0 max-gamma)] [1.0 max-h]]
        semi-variogram-fn (if (fn? semi-variogram)
                            (rbf->variogram semi-variogram)
                            (semi-variograms semi-variogram))
        target (fn [& r] (objective-least-squares
                         (semi-variogram-fn (zipmap [:nugget :sill :range :beta] r))
                         empirical-variogram))]
    [semi-variogram-fn (->> (optim/scan-and-minimize :lbfgsb target
                                                     {:bounds (if (#{:power :expower} semi-variogram)
                                                                (conj bounds [-8.0 8.0])
                                                                bounds)})
                            (first)
                            (zipmap [:nugget :sill :range :beta]))]))

(defn fit
  ([xss ys semi-variogram] (fit xss ys semi-variogram nil))
  ([xss ys semi-variogram opts] (fit (empirical xss ys opts) semi-variogram))
  ([empirical-variogram semi-variogram]
   (let [[f p] (fit-params empirical-variogram semi-variogram)]
     (f p))))
