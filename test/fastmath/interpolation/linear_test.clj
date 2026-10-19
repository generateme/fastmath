(ns fastmath.interpolation.linear-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.interpolation.linear :as sut]))

(set! *warn-on-reflection* true)

;; Oracles: numpy.interp (inside) and scipy interp1d with fill_value='extrapolate' (all queries),
;; see utils/fastmath/dev/generate_interp1d_reference.py.

(def ^:private ref-data
  (delay (edn/read-string (slurp (io/resource "interpolation/interp1d_reference.edn")))))

(defn- property
  "Returns the :property of the ex-info (:type :invalid-input) thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (= :invalid-input (:type d)) (:property d))))))

(defn- ymax ^double [ys] (double (reduce max 0.0 (map #(Math/abs (double %)) ys))))

(t/deftest linear-edge-cases
  (t/testing "two points: the straight line, also outside"
    (let [f (sut/linear [0 2] [1 5])]
      (t/is (m/delta-eq 1.0 (f 0) 1e-15))
      (t/is (m/delta-eq 5.0 (f 2) 1e-15))
      (t/is (m/delta-eq 3.0 (f 1) 1e-15))
      (t/is (m/delta-eq -1.0 (f -1) 1e-15))
      (t/is (m/delta-eq 9.0 (f 4) 1e-15))))
  (t/testing "knots return the data"
    (let [xs [-3.0 -1.0 0.5 2.0 10.0]
          ys [4.0 -2.0 7.5 0.25 -8.0]
          f (sut/linear xs ys)]
      (doseq [[x y] (map vector xs ys)]
        (t/is (m/delta-eq y (f x) 1e-14)))))
  (t/testing "sequence types"
    (let [f1 (sut/linear [0 1 2] [0 1 4])
          f2 (sut/linear '(0 1 2) '(0 1 4))
          f3 (sut/linear (double-array [0 1 2]) (double-array [0 1 4]))
          f4 (sut/linear (map inc [-1 0 1]) (range 0 5 2))]
      (t/is (== 0.5 (f1 0.5) (f2 0.5) (f3 0.5)))
      (t/is (== 2.5 (f1 1.5)))
      (t/is (== 3.0 (f4 1.5)))))
  (t/testing "NaN query gives NaN, NaN and infinity in ys are allowed and propagate to the neighbouring segments"
    (let [f (sut/linear [0 1 2 3] [0 ##NaN 2 3])]
      (t/is (m/nan? (f ##NaN)))
      (t/is (m/nan? (f 0)) "a node next to a NaN gives NaN (0 * NaN)")
      (t/is (m/nan? (f 0.5)))
      (t/is (m/nan? (f 1.5)))
      (t/is (m/delta-eq 2.5 (f 2.5) 1e-15))))
  (t/testing "a vector of vectors is not a valid y"
    (t/is (thrown? Exception (sut/linear [0 1] [[0] [1]])))))

(t/deftest linear-validation
  (t/testing "lengths"
    (t/is (= :length-mismatch (property #(sut/linear [0 1 2] [0 1]))))
    (t/is (= :length-mismatch (property #(sut/linear [0 1] [0 1 2])))))
  (t/testing "number of points"
    (t/is (= :too-few-points (property #(sut/linear [] []))))
    (t/is (= :too-few-points (property #(sut/linear [1] [1]))))
    (t/is (nil? (property #(sut/linear [1 2] [1 2])))))
  (t/testing "finite xs"
    (t/is (= :not-finite (property #(sut/linear [0 ##NaN 2] [0 1 2]))))
    (t/is (= :not-finite (property #(sut/linear [0 1 ##Inf] [0 1 2]))))
    (t/is (= :not-finite (property #(sut/linear [##-Inf 1 2] [0 1 2])))))
  (t/testing "strictly increasing xs: unsorted, duplicates, descending"
    (t/is (= :not-increasing (property #(sut/linear [0 2 1] [0 1 2]))))
    (t/is (= :not-increasing (property #(sut/linear [0 1 1] [0 1 2]))))
    (t/is (= :not-increasing (property #(sut/linear [1 1] [0 1]))))
    (t/is (= :not-increasing (property #(sut/linear [2 1 0] [0 1 2])))))
  (t/testing "ex-data names the method"
    (t/is (= :linear (:method (try (sut/linear [1] [1]) (catch clojure.lang.ExceptionInfo e (ex-data e))))))))

(t/deftest linear-reference
  ;; measured: 1.1e-16 * max(ys) inside, 5.5e-15 relative outside
  (doseq [d (:linear @ref-data)]
    (let [f (sut/linear (:xs d) (:ys d))
          lo (first (:xs d))
          hi (last (:xs d))
          sc (ymax (:ys d))]
      (doseq [[q e] (map vector (:qs d) (:expected d))]
        (t/is (m/delta-eq e (f q) (* 1e-12 sc) 1e-12) (str (:name d) " at " q (if (<= lo q hi) " inside" " outside")))))))
