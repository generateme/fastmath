(ns fastmath.interpolation.akima-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.interpolation.akima :as sut])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; Reference values: utils/fastmath/dev/generate_akima_loess_reference.py (scipy Akima1DInterpolator,
;; methods akima and makima, extrapolate=True).
(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "interpolation/akima_loess_reference.edn")))))

(defn- property
  "Returns the :property of the invalid-input ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :akima (:method d))) (:property d))))))

(defn- fv
  ^double [f x]
  (double (f x)))

(defn- max-err
  "Maximum absolute difference between `(f q)` and the expected values."
  ^double [f queries expected]
  (reduce (fn [^double acc [q e]]
            (m/max acc (m/abs (m/- (fv f q) (double e)))))
          0.0 (map vector queries expected)))

(defn- ymax
  ^double [ys]
  (reduce (fn [^double a y] (m/max a (m/abs (double y)))) 0.0 ys))

(defn- fourth-difference
  "Fourth difference of 5 samples of `f` with the step `h` starting at `a`: zero for a cubic polynomial."
  ^double [f ^double a ^double h]
  (let [v (double-array (map (fn [i] (fv f (m/+ a (m/* h (double i))))) (range 5)))]
    (m/+ (Array/aget v 0) (m/* -4.0 (Array/aget v 1)) (m/* 6.0 (Array/aget v 2))
         (m/* -4.0 (Array/aget v 3)) (Array/aget v 4))))

(def ^:private xs5 [0.0 1.0 2.0 3.5 4.0])
(def ^:private ys5 [2.0 -1.0 0.5 3.0 1.0])

(t/deftest akima-edge-cases
  (t/testing "too few points: 0, 1, 2, 3, 4 (scipy accepts 2), exactly the minimum 5 works"
    (t/is (= :too-few-points (property #(sut/akima [] []))))
    (t/is (= :too-few-points (property #(sut/akima [0.0] [1.0]))))
    (t/is (= :too-few-points (property #(sut/akima [0.0 1.0] [1.0 2.0]))))
    (t/is (= :too-few-points (property #(sut/akima [0.0 1.0 2.0] [1.0 2.0 3.0]))))
    (t/is (= :too-few-points (property #(sut/akima [0.0 1.0 2.0 3.0] [1.0 2.0 3.0 4.0]))))
    (t/is (= :too-few-points (property #(sut/akima [0.0 1.0 2.0 3.0] [1.0 2.0 3.0 4.0] {:method :makima}))))
    (t/is (fn? (sut/akima xs5 ys5))))
  (t/testing "invalid xs"
    (t/is (= :length-mismatch (property #(sut/akima xs5 [1.0 2.0 3.0 4.0]))))
    (t/is (= :length-mismatch (property #(sut/akima [0.0 1.0 2.0 3.0 4.0 5.0] ys5))))
    (t/is (= :not-finite (property #(sut/akima [0.0 1.0 ##NaN 3.0 4.0] ys5))))
    (t/is (= :not-finite (property #(sut/akima [0.0 1.0 2.0 3.0 ##Inf] ys5))))
    (t/is (= :not-finite (property #(sut/akima [##-Inf 1.0 2.0 3.0 4.0] ys5))))
    (t/is (= :not-increasing (property #(sut/akima [0.0 2.0 1.0 3.0 4.0] ys5))))
    (t/is (= :not-increasing (property #(sut/akima [4.0 3.0 2.0 1.0 0.0] ys5))))
    (t/is (= :not-increasing (property #(sut/akima [0.0 1.0 1.0 3.0 4.0] ys5))))
    (t/is (= :not-increasing (property #(sut/akima [0.0 1.0 2.0 3.0 3.0] ys5)))))
  (t/testing "unknown method"
    (t/is (= :unknown-method (property #(sut/akima xs5 ys5 {:method :pchip}))))
    (t/is (= :unknown-method (property #(sut/akima xs5 ys5 {:method "akima"})))))
  (t/testing "nil options and nil method mean the default"
    (let [f (sut/akima xs5 ys5)
          expected (fv f 1.7)]
      (t/is (= expected (fv (sut/akima xs5 ys5 nil) 1.7)))
      (t/is (= expected (fv (sut/akima xs5 ys5 {}) 1.7)))
      (t/is (= expected (fv (sut/akima xs5 ys5 {:method nil}) 1.7)))
      (t/is (= expected (fv (sut/akima xs5 ys5 {:method :akima}) 1.7)))))
  (t/testing "NaN ys are allowed at construction"
    (t/is (fn? (sut/akima xs5 [1.0 ##NaN 2.0 3.0 4.0])))))

(t/deftest akima-nodes-and-queries
  (doseq [method [:akima :makima]
          :let [f (sut/akima xs5 ys5 {:method method})]]
    (t/testing (str method " passes through the nodes")
      (doseq [[x y] (map vector xs5 ys5)]
        (t/is (m/delta-eq y (fv f x) 1.0e-13))))
    (t/testing (str method " NaN query")
      (t/is (m/nan? (fv f ##NaN))))
    (t/testing (str method " long and double queries, far outside the domain stay finite")
      (t/is (= (fv f 2) (fv f 2.0)))
      (t/is (m/valid-double? (fv f 1.0e12)))
      (t/is (m/valid-double? (fv f -1.0e12))))))

(t/deftest akima-linear-data
  ;; four equal consecutive slopes: both weights vanish, the mean fallback gives the exact line
  (let [xs [-2.0 -1.5 -0.3 0.2 1.1 2.0 3.7]
        line (fn ^double [^double x] (m/+ 1.0 (m/* 2.0 x)))]
    (doseq [method [:akima :makima]
            :let [f (sut/akima xs (map line xs) {:method method})]
            x [-10.0 -2.0 -1.0 0.0 0.2 1.5 3.7 8.0]]
      (t/is (m/delta-eq (line x) (fv f x) 1.0e-12)))))

(t/deftest akima-methods-differ
  (let [xs [0.0 0.7 1.1 2.0 2.4 3.9 4.0 5.5 6.1 7.3]
        ys [0.0 3.0 -2.0 4.0 -1.0 5.0 -3.0 2.0 0.0 1.0]
        fa (sut/akima xs ys)
        fm (sut/akima xs ys {:method :makima})]
    (t/is (> (double (reduce (fn [^double a ^double x] (m/max a (m/abs (m/- (fv fa x) (fv fm x)))))
                             0.0 (range 0.0 7.3 0.1)))
             1.0e-2))))

(t/deftest akima-extrapolation-continues-end-cubic
  ;; the end cubic is a cubic polynomial: 5 equidistant samples outside the domain have zero 4th difference
  (doseq [method [:akima :makima]
          :let [f (sut/akima xs5 ys5 {:method method})]]
    (t/is (m/delta-eq 0.0 (fourth-difference f -4.0 0.5) 1.0e-11))
    (t/is (m/delta-eq 0.0 (fourth-difference f 4.5 0.5) 1.0e-11))
    (t/testing "sanity: inside the domain the curve is not one cubic"
      (t/is (> (m/abs (fourth-difference f 0.0 1.0)) 1.0e-3)))))

(t/deftest akima-nan-ys
  ;; a NaN value changes the slopes of the nodes bad-2 ... bad+2 only; the NaN of the last node must not poison the rest
  (let [xs (mapv double (range 15))
        ys (mapv (fn [x] (m/sin (double x))) xs)]
    (doseq [[bad-idx far] [[7 [0.0 0.5 1.5 2.5 3.5 11.5 12.5 13.5 14.0 15.5 -1.0]]
                           [13 [0.0 0.5 1.5 2.5 5.5 8.5 9.5 -1.0]]
                           [1 [4.5 6.5 9.5 12.5 14.0 15.5]]]
            method [:akima :makima]
            :let [fb (sut/akima xs (assoc ys bad-idx ##NaN) {:method method})
                  fz (sut/akima xs (assoc ys bad-idx 0.0) {:method method})
                  b (double bad-idx)]]
      (t/testing (str method " NaN at " bad-idx " spreads only to the neighbourhood of the bad node")
        (t/is (m/nan? (fv fb b)))
        (t/is (m/nan? (fv fb (m/- b 0.5))))
        (t/is (m/nan? (fv fb (m/min 14.0 (m/+ b 0.5)))))
        (doseq [x far]
          (t/is (m/delta-eq (fv fz x) (fv fb x) 1.0e-12) (str "x=" x)))))))

(t/deftest akima-vs-scipy
  (doseq [{:keys [name xs ys queries] :as d} (:akima @reference)
          method [:akima :makima]
          :let [f (sut/akima xs ys {:method method})
                tol (* 1.0e-12 (ymax ys))]]
    ;; max abs error 1.4e-14 (extrema), relative to the largest |y| at most 1.2e-15 on all 12 datasets
    (t/is (<= (max-err f queries (get d method)) tol) (str name " " method))))

(t/deftest akima-end-intervals-and-extrapolation-vs-scipy
  ;; the old ACM implementation differed from scipy by up to 0.25 here
  (doseq [{:keys [name xs ys queries] :as d} (:akima @reference)
          method [:akima :makima]
          :let [f (sut/akima xs ys {:method method})
                tol (* 1.0e-12 (ymax ys))
                n (count xs)
                lo (double (first xs))
                hi (double (last xs))
                ;; the first 40 queries are a grid reaching 25 percent beyond both ends
                expected (zipmap queries (get d method))
                end-q (filter (fn [^double q] (or (m/< q (double (second xs))) (m/> q (double (nth xs (- n 2))))))
                              queries)
                outside (filter (fn [^double q] (or (m/< q lo) (m/> q hi))) (take 40 queries))]]
    (t/is (seq outside) name)
    (t/is (seq end-q) name)
    (doseq [q (concat end-q outside)]
      (t/is (m/delta-eq (double (get expected q)) (fv f q) tol) (str name " " method " x=" q)))))

(t/deftest akima-fewer-points-than-required
  ;; Suggestion (not enforced, the Spec requires 5 points): the algorithm itself works for 2, 3 and 4 points
  ;; like scipy (2 points: linear). Measured max error against scipy 1.8e-15.
  (doseq [{:keys [name xs ys queries] :as d} (:akima-small @reference)
          method [:akima :makima]
          :let [f (#'sut/hermite (double-array xs) (double-array ys) (= method :makima))]]
    (t/is (<= (max-err f queries (get d method)) 1.0e-13) (str name " " method))))
