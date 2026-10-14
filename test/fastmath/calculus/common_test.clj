(ns fastmath.calculus.common-test
  (:require [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.calculus.common :as sut])
  (:import [fastmath.calculus.common Estimate]))

(deftype Region [^double v ^double e]
  Estimate
  (value [_] v)
  (error [_] e))

(defn- drain
  "Errors of all queued regions in the order of polling."
  [^java.util.PriorityQueue q]
  (loop [res []]
    (if-let [^Estimate x (.poll q)]
      (recur (conj res (.error x)))
      res)))

(t/deftest max-error-queue
  (t/testing "empty queue"
    (t/is (nil? (.poll (sut/max-error-queue))))
    (t/is (= [0.0 0.0] (sut/estimates-sum (sut/max-error-queue)))))
  (t/testing "one region"
    (let [q (sut/max-error-queue)]
      (.add q (Region. 2.0 0.5))
      (t/is (= [0.5] (drain q)))))
  (t/testing "polling returns non-increasing errors for many random regions"
    (let [rng (r/rng :jdk 1)
          q (sut/max-error-queue)
          errors (repeatedly 10000 #(r/drandom rng))]
      (doseq [e errors] (.add q (Region. 1.0 e)))
      (t/is (= (reverse (sort errors)) (drain q)))))
  (t/testing "regions with equal errors are all kept (polynomials give many zero errors)"
    (let [q (sut/max-error-queue)]
      (doseq [e [0.0 0.0 0.0 1.0 1.0]] (.add q (Region. 1.0 e)))
      (t/is (= 5 (.size q)))
      (t/is (= [1.0 1.0 0.0 0.0 0.0] (drain q)))))
  (t/testing "a NaN error is the largest, an infinite one comes after it"
    (let [q (sut/max-error-queue)]
      (doseq [e [1.0 ##Inf ##NaN 3.0 ##-Inf]] (.add q (Region. 1.0 e)))
      (let [res (drain q)]
        (t/is (m/nan? (first res)))
        (t/is (= [##Inf 3.0 1.0 ##-Inf] (rest res))))))
  (t/testing "sums of values and errors"
    (let [q (sut/max-error-queue)]
      (doseq [[v e] [[1.0 0.1] [2.0 0.2] [-0.5 0.3]]] (.add q (Region. v e)))
      (let [[v e] (sut/estimates-sum q)]
        (t/is (m/delta-eq 2.5 v 1.0e-15))
        (t/is (m/delta-eq 0.6 e 1.0e-15))))))

(t/deftest normalize-bounds
  (t/are [lower upper res] (= res (sut/normalize-bounds lower upper))
    0.0 1.0 [1.0 0.0 1.0]
    1.0 0.0 [-1.0 0.0 1.0]
    2.0 2.0 [1.0 2.0 2.0]
    -1.0 -2.0 [-1.0 -2.0 -1.0]
    ##Inf 0.0 [-1.0 0.0 ##Inf]
    0.0 ##-Inf [-1.0 ##-Inf 0.0]
    ##Inf ##-Inf [-1.0 ##-Inf ##Inf]
    ##-Inf ##Inf [1.0 ##-Inf ##Inf])
  (t/testing "NaN bounds and two equal infinite bounds throw ex-info with the bounds"
    (doseq [[lower upper] [[##NaN 0.0] [0.0 ##NaN] [##NaN ##NaN] [##Inf ##Inf] [##-Inf ##-Inf]]]
      (let [e (try (sut/normalize-bounds lower upper) nil (catch clojure.lang.ExceptionInfo e e))]
        (t/is (some? e))
        (t/is (contains? (ex-data e) :lower))
        (t/is (contains? (ex-data e) :upper))))))

(t/deftest normalize-bounds-multi
  (t/are [lower upper res] (= res (sut/normalize-bounds-multi lower upper))
    [0 0] [1 1] [1.0 [0.0 0.0] [1.0 1.0]]
    [1 0] [0 1] [-1.0 [0.0 0.0] [1.0 1.0]]
    [1 1] [0 0] [1.0 [0.0 0.0] [1.0 1.0]]
    [1 1 1] [0 0 0] [-1.0 [0.0 0.0 0.0] [1.0 1.0 1.0]]
    [##Inf 0] [0 1] [-1.0 [0.0 0.0] [##Inf 1.0]]
    [2] [2] [1.0 [2.0] [2.0]])
  (t/testing "empty and mismatched bounds throw ex-info"
    (doseq [[lower upper] [[[] []] [[1 0] [0]] [[1] [0 1]] [nil nil]]]
      (t/is (thrown? clojure.lang.ExceptionInfo (sut/normalize-bounds-multi lower upper)))))
  (t/testing "a NaN bound and two equal infinite bounds throw ex-info in any dimension"
    (doseq [[lower upper] [[[##NaN 0] [1 1]] [[0 0] [1 ##NaN]] [[0 ##Inf] [1 ##Inf]] [[0 ##-Inf] [1 ##-Inf]]]]
      (t/is (thrown? clojure.lang.ExceptionInfo (sut/normalize-bounds-multi lower upper))))))

(t/deftest drop-nil-options
  (t/is (= {} (sut/drop-nil-options nil)))
  (t/is (= {} (sut/drop-nil-options {})))
  (t/is (= {:a false :b 0 :c ##Inf} (sut/drop-nil-options {:a false :b 0 :c ##Inf :d nil :e nil})) "false is not nil"))

(t/deftest orient
  (t/is (= 3.0 (sut/orient 1.0 3.0)))
  (t/is (= -3.0 (sut/orient -1.0 3.0)))
  (t/is (= {:result -3.0 :error 0.5} (sut/orient -1.0 {:result 3.0 :error 0.5})))
  (t/is (= {:result 3.0 :error 0.5} (sut/orient 1.0 {:result 3.0 :error 0.5})))
  (t/is (m/nan? (sut/orient -1.0 ##NaN))))

(t/deftest subst-1d-inverse
  ;; measured: x -> t -> x within 5.1e-14 (relative to max(1, |x|)) for |x| <= 1e3, t -> x -> t within 4.4e-16
  (t/are [lower upper xs] (every? (fn [^double x]
                                    (let [[_ _ _] (sut/subst-1d identity lower upper)
                                          t ((sut/subst-1d-inverse lower upper) x)
                                          x' (cond
                                               (and (m/neg-inf? lower) (m/pos-inf? upper)) (m// t (m/- 1.0 (m/* t t)))
                                               (m/neg-inf? lower) (m/+ upper (m// t (m/inc t)))
                                               (m/pos-inf? upper) (m/+ lower (m// t (m/- 1.0 t)))
                                               :else t)]
                                      (m/< (m/abs (m/- x x')) (m/* 1.0e-12 (m/max 1.0 (m/abs x))))))
                                  xs)
    ##-Inf ##Inf [-1000.0 -3.0 -1.0e-3 0.0 0.5 7.0 1000.0]
    ##-Inf 2.0 [-1000.0 -3.0 0.0 2.0]
    1.5 ##Inf [1.5 2.0 100.0 1000.0]
    0.0 1.0 [0.0 0.3 1.0])
  (t/testing "the image of the interval is the interval of the substitution"
    (let [lower 1.5
          t ((sut/subst-1d-inverse lower ##Inf) lower)]
      (t/is (= 0.0 t)))
    (let [t ((sut/subst-1d-inverse ##-Inf 2.0) 2.0)]
      (t/is (= 0.0 t)))))
