(ns fastmath.interpolation.bootstrap-test
  (:require [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.interpolation :as interp]
            [fastmath.interpolation.bootstrap :as sut]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- property
  "The :property of the invalid-input ex-info thrown by `f`, nil when it doesn't throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)] (when (= :invalid-input (:type d)) (:property d))))))

(def ^:private rng-data (r/rng :mersenne 11))
(def ^:private xs (vec (range 0.0 20.0)))
(def ^:private ys (mapv (fn [^double x] (m/+ (m/* 2.0 x) (r/grandom rng-data 0.0 1.0))) xs))
(def ^:private line (fn [xs ys] (let [{:keys [beta intercept]} (fastmath.ml.regression/lm ys (map vector xs))]
                                  (fn [^double x] (m/+ (double intercept) (m/* (double (first beta)) x))))))

(t/deftest bootstrap-band
  (let [b (sut/bootstrap line xs ys {:samples 300 :rng (r/rng :mersenne 1)})]
    (t/testing "record, fn?, the fit of the data"
      (t/is (fn? b))
      (t/is (instance? fastmath.interpolation.bootstrap.Bootstrap b))
      (t/is (m/delta-eq (double ((line xs ys) 7.5)) (double (b 7.5)))))
    (t/testing "the band contains the true line at most points and narrows with the level"
      (t/is (>= (count (filter (fn [^double x] (let [[lo hi] (b x 0.95)] (m/<= (double lo) (m/* 2.0 x) (double hi)))) xs)) 17))
      (let [[lo95 hi95] (b 10.0 0.95) [lo50 hi50] (b 10.0 0.5)]
        (t/is (m/< (double lo95) (double lo50) (double hi50) (double hi95)))))
    (t/testing "every :ci gives an interval around the fit; :estimation is passed"
      (doseq [ci [:percentile :basic :normal :bc]]
        (let [bc (sut/bootstrap line xs ys {:samples 300 :rng (r/rng :mersenne 1) :ci ci})
              [lo hi] (bc 10.0 0.9)]
          (t/is (m/< (double lo) (double (bc 10.0)) (double hi)) (str ci))))
      (let [bp (sut/bootstrap line xs ys {:samples 300 :rng (r/rng :mersenne 1) :estimation :r7})]
        (t/is (not= (bp 10.0 0.9) (b 10.0 0.9)) "a different quantile strategy")))
    (t/testing "a fixed seed reproduces the band"
      (t/is (= (b 3.0 0.9) ((sut/bootstrap line xs ys {:samples 300 :rng (r/rng :mersenne 1)}) 3.0 0.9))))))

(t/deftest bootstrap-data-and-failures
  (t/testing "unsorted data sorted, repeated x averaged"
    (let [b (sut/bootstrap (fn [xs ys] (let [m (zipmap xs ys)] (fn [x] (double (m (double x) ##NaN)))))
                           [2.0 1.0 2.0] [5.0 1.0 3.0] {:samples 1})]
      (t/is (= 1.0 (b 1.0)))
      (t/is (= 4.0 (b 2.0)) "mean of 5 and 3")))
  (t/testing "failing refits are skipped and counted; none left - [fit fit]"
    (let [xs6 [0.0 1.0 2.0 3.0 4.0 5.0]
          b (interp/bootstrap :akima xs6 [0.0 1.0 0.0 1.0 0.0 1.0] {:samples 100 :rng (r/rng :mersenne 3)})]
      (t/is (m/pos? (long (:failed b))) "akima needs 5 distinct points")
      (t/is (= (+ (count (:fits b)) (long (:failed b))) 100)))
    (let [calls (atom 0)
          ;; the first call fits the data, every refit throws
          b (sut/bootstrap (fn [_ _] (if (= 1 (swap! calls inc)) (fn [_] 4.0) (throw (ex-info "no refit" {}))))
                           [0.0 1.0 2.0] [1.0 2.0 3.0] {:samples 5})]
      (t/is (= 5 (:failed b)))
      (t/is (= [4.0 4.0] (b 1.0 0.9)))))
  (t/testing "one sample, one point"
    (let [b (sut/bootstrap (fn [_ ys] (fn [_] (double (first ys)))) [1.0] [7.0] {:samples 1})]
      (t/is (= [7.0 7.0] (b 0.0 0.5))))))

(t/deftest bootstrap-interpolation-wrapper
  (t/testing "method keywords with params, functions"
    (let [b (interp/bootstrap :loess xs ys {:samples 20 :rng (r/rng :mersenne 5) :params [{:bandwidth 0.5}]})]
      (t/is (m/delta-eq (double ((interp/loess xs ys {:bandwidth 0.5}) 4.2)) (double (b 4.2))))
      (t/is (= 2 (count (b 4.2 0.9)))))
    (t/is (fn? (interp/bootstrap line xs ys {:samples 2}))))
  (t/testing "invalid input"
    (t/is (= :method (property #(interp/bootstrap :no-such-method xs ys))))
    (t/is (= :method (property #(interp/bootstrap 42 xs ys))))
    (t/is (= :samples (property #(sut/bootstrap line xs ys {:samples 0}))))
    (t/is (= :samples (property #(sut/bootstrap line xs ys {:samples 2.5}))))
    (t/is (= :ci (property #(sut/bootstrap line xs ys {:ci :bca}))))
    (t/is (= :fit-fn (property #(sut/bootstrap (fn [_ _] 42) xs ys {:samples 2}))))
    (t/is (= :length-mismatch (property #(sut/bootstrap line xs (rest ys)))))
    (let [b (sut/bootstrap line xs ys {:samples 2})]
      (t/is (= :level (property #(b 1.0 0.0))))
      (t/is (= :level (property #(b 1.0 1.0))))
      (t/is (= :level (property #(b 1.0 -0.5)))))))
