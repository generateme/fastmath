(ns fastmath.interpolation.loess-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.interpolation.loess :as sut]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; Reference values: utils/fastmath/dev/generate_akima_loess_reference.py
;; statsmodels.nonparametric.lowess(delta=0), the natural cubic spline of scipy CubicSpline through the smoothed values
;; and, for user weights, the Cython routine _smoothers_lowess.lowess with resid_weights (iters 0).
(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "interpolation/akima_loess_reference.edn")))))

(def ^:private datasets
  (delay (into {} (map (juxt :name identity)) (:loess-data @reference))))

(defn- property
  "Returns the :property of the invalid-input ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :loess (:method d))) (:property d))))))

(defn- fv
  ^double [f x]
  (double (f x)))

(defn- max-err
  "Maximum absolute difference between `(f q)` and the expected values."
  ^double [f queries expected]
  (reduce (fn [^double acc [q e]]
            (m/max acc (m/abs (m/- (fv f q) (double e)))))
          0.0 (map vector queries expected)))

(def ^:private xs7 [0.0 0.5 1.5 2.0 3.5 4.0 6.0])
(def ^:private ys7 [1.0 2.0 0.0 3.0 2.5 4.0 1.0])

(t/deftest loess-window-rule
  (t/testing "k = int(bandwidth * n + 1e-10), at least 2"
    (t/is (= 2 (#'sut/window-size 0.3 7)))
    (t/is (= 1 (#'sut/window-size 0.3 6)))
    (t/is (= 15 (#'sut/window-size 0.3 50)))
    (t/is (= 4 (#'sut/window-size 0.3 15)))
    (t/is (= 7 (#'sut/window-size 1.0 7)))
    (t/is (= 2 (#'sut/window-size (/ 2.0 7.0) 7)) "rounding up when the product is just below an integer")
    (t/is (= 1 (#'sut/window-size 0.28 7)))))

(t/deftest loess-edge-cases
  (t/testing "points: empty, one, two"
    (t/is (= :too-few-points (property #(sut/loess [] [] {:bandwidth 1.0}))))
    (t/is (= :too-few-points (property #(sut/loess [1.0] [2.0] {:bandwidth 1.0}))))
    (t/is (= :bandwidth (property #(sut/loess [0.0 1.0] [2.0 3.0]))) "default bandwidth gives k = 0 for two points")
    (let [f (sut/loess [0.0 1.0] [2.0 3.0] {:bandwidth 1.0})]
      (t/is (m/delta-eq 2.0 (fv f 0.0) 1.0e-13))
      (t/is (m/delta-eq 3.0 (fv f 1.0) 1.0e-13))
      (t/is (m/delta-eq 4.0 (fv f 2.0) 1.0e-12) "linear continuation"))
    (t/is (fn? (sut/loess [0.0 1.0 2.0] [2.0 3.0 1.0] {:bandwidth 1.0}))))
  (t/testing "minimal count for the default bandwidth 0.3 is 7 points"
    (t/is (= :bandwidth (property #(sut/loess (range 6.0) (range 6.0)))))
    (t/is (fn? (sut/loess (range 7.0) (range 7.0)))))
  (t/testing "bandwidth range and the minimal window"
    (t/is (= :bandwidth (property #(sut/loess xs7 ys7 {:bandwidth 0.0}))))
    (t/is (= :bandwidth (property #(sut/loess xs7 ys7 {:bandwidth -0.5}))))
    (t/is (= :bandwidth (property #(sut/loess xs7 ys7 {:bandwidth 1.5}))))
    (t/is (= :bandwidth (property #(sut/loess xs7 ys7 {:bandwidth ##NaN}))))
    (t/is (= :bandwidth (property #(sut/loess xs7 ys7 {:bandwidth 0.28}))))
    (t/is (= :bandwidth (property #(sut/loess xs7 ys7 {:bandwidth "0.3"}))))
    (t/is (fn? (sut/loess xs7 ys7 {:bandwidth (/ 2.0 7.0)})))
    (t/is (fn? (sut/loess xs7 ys7 {:bandwidth 1.0}))))
  (t/testing "xs"
    (t/is (= :length-mismatch (property #(sut/loess xs7 [1.0 2.0]))))
    (t/is (= :not-increasing (property #(sut/loess [0.0 2.0 1.0 3.0 4.0 5.0 6.0] ys7))))
    (t/is (= :not-increasing (property #(sut/loess [0.0 1.0 1.0 3.0 4.0 5.0 6.0] ys7))))
    (t/is (= :not-finite (property #(sut/loess [0.0 1.0 ##NaN 3.0 4.0 5.0 6.0] ys7))))
    (t/is (= :not-finite (property #(sut/loess [0.0 1.0 2.0 3.0 4.0 5.0 ##Inf] ys7)))))
  (t/testing "iters, accuracy, weights"
    (t/is (= :iters (property #(sut/loess xs7 ys7 {:bandwidth 0.5 :iters -1}))))
    (t/is (= :accuracy (property #(sut/loess xs7 ys7 {:bandwidth 0.5 :accuracy -1.0e-3}))))
    (t/is (= :length-mismatch (property #(sut/loess xs7 ys7 {:bandwidth 0.5 :weights [1.0 1.0]}))))
    (t/is (= :invalid-weights (property #(sut/loess xs7 ys7 {:bandwidth 0.5 :weights [1.0 1.0 -1.0 1.0 1.0 1.0 1.0]}))))
    (t/is (= :invalid-weights (property #(sut/loess xs7 ys7 {:bandwidth 0.5 :weights [1.0 1.0 ##NaN 1.0 1.0 1.0 1.0]}))))
    (t/is (fn? (sut/loess xs7 ys7 {:bandwidth 0.5 :iters 0 :accuracy 0.0 :weights [1.0 0.0 1.0 1.0 1.0 1.0 1.0]}))))
  (t/testing "nil options mean the defaults"
    (let [xs (range 20.0)
          ys (map (fn [^double x] (m/sin x)) xs)
          expected (fv (sut/loess xs ys) 7.3)]
      (t/is (= expected (fv (sut/loess xs ys nil) 7.3)))
      (t/is (= expected (fv (sut/loess xs ys {:bandwidth nil :iters nil :accuracy nil :weights nil}) 7.3)))
      (t/is (= expected (fv (sut/loess xs ys {:bandwidth 0.3 :iters 2 :accuracy 1.0e-12}) 7.3)))))
  (t/testing "NaN query gives NaN, far outside the domain stays finite"
    (let [f (sut/loess xs7 ys7 {:bandwidth 0.5})]
      (t/is (m/nan? (fv f ##NaN)))
      (t/is (m/valid-double? (fv f 1.0e9)))
      (t/is (m/valid-double? (fv f -1.0e9))))))

(t/deftest loess-exact-data
  (t/testing "k = 2 keeps the data (the farthest neighbour has the zero weight)"
    (let [f (sut/loess xs7 ys7)
          f2 (sut/loess xs7 ys7 {:bandwidth 0.3 :iters 0})]
      (doseq [[x y] (map vector xs7 ys7)]
        (t/is (m/delta-eq y (fv f x) 1.0e-12))
        (t/is (m/delta-eq y (fv f2 x) 1.0e-12)))))
  (t/testing "constant and linear data are reproduced for every iteration count"
    (let [xs (mapv (fn [^double i] (m/+ i (m/* 0.3 (m/sin (m/* 2.0 i))))) (range 25.0))]
      (doseq [iters [0 1 2 5]
              bw [0.2 0.5 1.0]
              :let [fc (sut/loess xs (repeat 25 3.5) {:bandwidth bw :iters iters})
                    fl (sut/loess xs (map (fn [^double x] (m/+ 1.0 (m/* -2.0 x))) xs) {:bandwidth bw :iters iters})]
              ^double x xs]
        ;; local linear regression is exact for lines; observed error below 1e-13
        (t/is (m/delta-eq 3.5 (fv fc x) 1.0e-12))
        (t/is (m/delta-eq (m/+ 1.0 (m/* -2.0 x)) (fv fl x) 1.0e-12))))))

(t/deftest loess-structural-properties
  (let [{:keys [xs ys]} (@datasets "nonuniform50")]
    (t/testing "scaling of x and y does not change the smoothed values (relative)"
      (let [f (sut/loess xs ys {:bandwidth 0.3})]
        (doseq [[sx sy] [[1.0e-9 1.0] [1.0e9 1.0] [1.0 1.0e6] [1.0 1.0e-6] [1.0e-4 1.0e4]]
                :let [fs (sut/loess (map (fn [^double x] (m/* x (double sx))) xs)
                                    (map (fn [^double y] (m/* y (double sy))) ys)
                                    {:bandwidth 0.3})]
                x xs
                :let [e (m/* (fv f x) (double sy))]]
          ;; observed relative error 1e-14 (the lower bound 1e-12 of the regression variance is relative to the radius)
          (t/is (m/delta-eq e (fv fs (m/* (double x) (double sx))) (m/* 1.0e-12 (m/max 1.0 (m/abs e)))) (str sx " " sy)))))
    (t/testing "a constant factor of all user weights changes nothing"
      (let [w (repeat 50 5.0)
            a (sut/loess xs ys {:bandwidth 0.3})
            b (sut/loess xs ys {:bandwidth 0.3 :weights w})]
        (doseq [x xs]
          (t/is (m/delta-eq (fv a x) (fv b x) 1.0e-12)))))
    (t/testing "accuracy stops the robustness iterations: a huge accuracy equals iters 0"
      (let [a (sut/loess xs ys {:bandwidth 0.3 :iters 0})
            b (sut/loess xs ys {:bandwidth 0.3 :iters 2 :accuracy 1.0e9})]
        (doseq [x xs]
          (t/is (= (fv a x) (fv b x))))))
    (t/testing "more robustness iterations move the fit away from the outliers"
      (let [rows (filter #(and (= "nonuniform50" (:name %)) (= 0.3 (:bandwidth %))) (:loess @reference))
            s0 (:smooth (first (filter #(= 0 (:iters %)) rows)))
            s2 (:smooth (first (filter #(= 2 (:iters %)) rows)))]
        ;; the outliers are y[5] (+8) and y[30] (-6)
        (t/is (> (m/abs (m/- (double (nth ys 5)) (double (nth s2 5))))
                 (m/abs (m/- (double (nth ys 5)) (double (nth s0 5))))))
        (t/is (> (m/abs (m/- (double (nth ys 30)) (double (nth s2 30))))
                 (m/abs (m/- (double (nth ys 30)) (double (nth s0 30))))))))))

(t/deftest loess-many-points
  (let [n 2000
        xs (mapv (fn [^double i] (m/+ i (m/* 0.4 (m/sin (m/* 3.0 i))))) (range (double n)))
        ;; k = 10 neighbours, the local linear bias is about 1e-4 (observed 1.2e-4)
        f (sut/loess xs (map (fn [^double x] (m/+ (m/sin (m/* 0.01 x)) 2.0)) xs) {:bandwidth 0.005})]
    (t/is (m/delta-eq (m/+ (m/sin 5.0) 2.0) (fv f 500.0) 1.0e-3))
    (t/is (m/valid-double? (fv f 2500.0)))))

(t/deftest loess-vs-statsmodels
  ;; smoothed values at the knots; max abs error 8.2e-15 over all 38 configurations (old ACM version: up to 1.45)
  (doseq [{:keys [name bandwidth iters smooth]} (:loess @reference)
          :let [{:keys [xs ys]} (@datasets name)]]
    (t/testing (str name " bandwidth " bandwidth " iters " iters)
      ;; accuracy 0.0 follows statsmodels exactly (also when the median residual is 0, k = 2)
      (t/is (<= (max-err (sut/loess xs ys {:bandwidth bandwidth :iters iters :accuracy 0.0}) xs smooth) 1.0e-12))
      ;; the default accuracy 1e-12 gives the same values on these noisy data
      (t/is (<= (max-err (sut/loess xs ys {:bandwidth bandwidth :iters iters}) xs smooth) 1.0e-12)))))

(t/deftest loess-spline-vs-scipy
  ;; natural cubic spline through the smoothed values, queries 10 percent beyond both ends:
  ;; max abs error 2.2e-12 outside the domain (end cubics), 4.7e-14 inside
  (doseq [{:keys [name bandwidth iters queries spline]} (:loess @reference)
          :let [{:keys [xs ys]} (@datasets name)
                f (sut/loess xs ys {:bandwidth bandwidth :iters iters})
                lo (double (first xs))
                hi (double (last xs))
                inside (filter (fn [[^double q _]] (and (m/>= q lo) (m/<= q hi))) (map vector queries spline))]]
    (t/is (<= (max-err f queries spline) 1.0e-11) (str name " " bandwidth " " iters))
    (t/is (<= (max-err f (map first inside) (map second inside)) 1.0e-12) (str name " " bandwidth " " iters))))

(t/deftest loess-user-weights-vs-statsmodels
  ;; Cython lowess with resid_weights = user weights (iters 0), one weight is exactly 0; max abs error 3.6e-15
  (doseq [{:keys [name bandwidth weights smooth]} (:loess-weighted @reference)
          :let [{:keys [xs ys]} (@datasets name)
                f (sut/loess xs ys {:bandwidth bandwidth :iters 0 :weights weights})]]
    (t/is (<= (max-err f xs smooth) 1.0e-12) (str name " " bandwidth))))
