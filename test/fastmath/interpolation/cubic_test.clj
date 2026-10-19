(ns fastmath.interpolation.cubic-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.interpolation.cubic :as sut]))

(set! *warn-on-reflection* true)

;; Oracle: scipy CubicSpline(bc_type='natural', extrapolate=True),
;; see utils/fastmath/dev/generate_interp1d_reference.py.

(def ^:private ref-data
  (delay (edn/read-string (slurp (io/resource "interpolation/interp1d_reference.edn")))))

(defn- property
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (= :invalid-input (:type d)) (:property d))))))

(defn- ymax ^double [ys] (double (reduce max 0.0 (map #(Math/abs (double %)) ys))))

(t/deftest cubic-edge-cases
  (t/testing "two points: the straight line, also outside"
    (let [f (sut/cubic [0 2] [1 5])]
      (t/is (m/delta-eq 1.0 (f 0) 1e-15))
      (t/is (m/delta-eq 5.0 (f 2) 1e-15))
      (t/is (m/delta-eq 3.0 (f 1) 1e-15))
      (t/is (m/delta-eq -3.0 (f -2) 1e-14))
      (t/is (m/delta-eq 9.0 (f 4) 1e-14))))
  (t/testing "three points: natural spline through the data, second derivative zero at the ends"
    ;; xs 0 1 2, ys 0 1 0: the natural spline is a cubic on each interval with y'' = 0 at 0 and 2,
    ;; y''(1) = -3: y = 1.5 x - 0.5 x^3 on [0, 1]
    (let [f (sut/cubic [0 1 2] [0 1 0])]
      (t/is (m/delta-eq 0.0 (f 0) 1e-15))
      (t/is (m/delta-eq 1.0 (f 1) 1e-15))
      (t/is (m/delta-eq 0.0 (f 2) 1e-15))
      (t/is (m/delta-eq (- (* 1.5 0.5) (* 0.5 0.125)) (f 0.5) 1e-14))
      (t/is (m/delta-eq (f 1.5) (f 0.5) 1e-14))))
  (t/testing "linear data is reproduced exactly, also in extrapolation"
    (let [xs [-2.0 -0.5 0.0 1.5 4.0]
          f (sut/cubic xs (map #(- (* 3.0 %) 1.0) xs))]
      (doseq [x [-10.0 -2.0 -1.0 0.3 2.0 4.0 7.0]]
        (t/is (m/delta-eq (- (* 3.0 x) 1.0) (f x) 1e-12)))))
  (t/testing "extrapolation continues the end cubic (it is not linear)"
    (let [f (sut/cubic [0 1 2 3] [0 1 0 2])
          ;; the cubic through four points of the first (last) interval, evaluated by Lagrange
          lagrange (fn [pts x]
                     (reduce + (for [xi pts]
                                 (* (f xi) (reduce * (for [xj pts :when (not= xi xj)]
                                                       (/ (- x xj) (- xi xj))))))))]
      ;; tolerance 1e-12: the Lagrange form amplifies the rounding errors
      (t/is (m/delta-eq (lagrange [0.0 0.2 0.5 0.9] -2.0) (f -2.0) 1e-12))
      (t/is (m/delta-eq (lagrange [2.1 2.4 2.7 3.0] 5.0) (f 5.0) 1e-12))
      ;; a linear continuation of the end slope would be different
      (t/is (not (m/delta-eq (+ (f 0.0) (* -2.0 (/ (- (f 1e-6) (f 0.0)) 1e-6))) (f -2.0) 1e-3)))))
  (t/testing "knots return the data"
    (let [xs [-3.0 -1.0 0.5 2.0 10.0]
          ys [4.0 -2.0 7.5 0.25 -8.0]
          f (sut/cubic xs ys)]
      (doseq [[x y] (map vector xs ys)]
        (t/is (m/delta-eq y (f x) 1e-13)))))
  (t/testing "NaN query gives NaN"
    (t/is (m/nan? ((sut/cubic [0 1 2] [0 1 4]) ##NaN))))
  (t/testing "sequence types"
    (t/is (m/delta-eq ((sut/cubic [0 1 2 3] [0 1 4 9]) 1.5)
                      ((sut/cubic (double-array [0 1 2 3]) (list 0 1 4 9)) 1.5) 1e-15))))

(t/deftest cubic-validation
  (t/is (= :length-mismatch (property #(sut/cubic [0 1 2] [0 1]))))
  (t/is (= :too-few-points (property #(sut/cubic [] []))))
  (t/is (= :too-few-points (property #(sut/cubic [1] [1]))))
  (t/is (= :not-finite (property #(sut/cubic [0 ##NaN 2] [0 1 2]))))
  (t/is (= :not-finite (property #(sut/cubic [0 1 ##Inf] [0 1 2]))))
  (t/testing "unsorted xs gave a silently wrong value, duplicates gave NaN"
    (t/is (= :not-increasing (property #(sut/cubic [0 2 1 3] [0 1 2 3]))))
    (t/is (= :not-increasing (property #(sut/cubic [0 1 1 3] [0 1 2 3]))))
    (t/is (= :not-increasing (property #(sut/cubic [3 2 1 0] [0 1 2 3])))))
  (t/is (nil? (property #(sut/cubic [0 1] [0 1])))))

(t/deftest cubic-reference
  ;; measured: 1.8e-15 * max(ys) inside, 7.1e-14 relative outside
  (doseq [d (:cubic @ref-data)]
    (let [f (sut/cubic (:xs d) (:ys d))
          lo (first (:xs d))
          hi (last (:xs d))
          sc (ymax (:ys d))]
      (doseq [[q e] (map vector (:qs d) (:expected d))]
        (t/is (m/delta-eq e (f q) (* 1e-12 sc) 1e-12) (str (:name d) " at " q (if (<= lo q hi) " inside" " outside")))))))

(t/deftest cubic-docstring
  (t/is (re-find #"natural" (:doc (meta #'sut/cubic))))
  (t/is (re-find #"(?i)extrapolation" (:doc (meta #'sut/cubic))))
  (t/is (re-find #"not-a-knot" (:doc (meta #'sut/cubic-2d)))))
